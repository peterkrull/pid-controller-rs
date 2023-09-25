#![feature(type_changing_struct_update)]
#![no_std]

use core::ops::{AddAssign, SubAssign};

use num_traits::float::Float;

/// The controller does not wrap
pub struct NoWrapping; 
/// The controller will wrap around from `min -> max` or `max -> min` gracefully
pub struct Wrapping<F> {
    min: F,
    max: F,
}

/// The controller has no output limit
pub struct NoOutputLimit;
/// The controller limits its output between `min` and `max`
pub struct OutputLimit<F> {
    min: F,
    max: F,
}

/// The controller does not apply anti-windup to integrator
pub struct NoAntiWindup;
/// The controller applies an anti windup method to the integrator
pub enum AntiWindup<F> {
    Conditional((F, F)),
}

/// The error is now low-pass filtered for the derivative controller
pub struct NoLpFilter;
/// The error is low-pass filtered with time constant `tau` for the derivative controller
pub struct LpFilter<F> {
    tau: F,
    signal: F,
}

/// Pid controller implementation
pub struct Pid<F, WRAP, LIMIT, WINDUP, LPFILT> {
    kp: F,
    ki: F,
    kd: F,
    ts: F,
    integral: F,
    integral_en: bool,
    prev_error: F,
    wrapping: WRAP,
    output_limit: LIMIT,
    anti_windup: WINDUP,
    lp_filter: LPFILT,
}

impl<F: Float + Default> Pid<F, NoWrapping, NoOutputLimit, NoAntiWindup, NoLpFilter> {
    /// ## PID controller
    /// 
    /// Construct a new PID controller
    /// 
    /// - `kp`: Proportional gain
    /// - `ki`: Integral gain
    /// - `kd`: Derivative gain
    /// - `ideal`: Controller on form `kp * ( 1 + ki/s + kd*s )`
    /// - `ts`: Sample time the controller is expected to run at
    pub fn new(
        kp: F,
        ki: F,
        kd: F,
        ideal: bool,
        ts: F,
    ) -> Pid<F, NoWrapping, NoOutputLimit, NoAntiWindup, NoLpFilter> {
        Pid {
            kp,
            ki: if ideal { ki * kp } else { ki },
            kd: if ideal { kd * kp } else { kd },
            ts,
            integral: F::zero(),
            integral_en: true,
            prev_error: F::zero(),
            wrapping: NoWrapping,
            output_limit: NoOutputLimit,
            anti_windup: NoAntiWindup,
            lp_filter: NoLpFilter,
        }
    }
}

impl<F: Float + Default, ANYWRAP, ANYLIMIT, ANYWINDUP, ANYFILT>
    Pid<F, ANYWRAP, ANYLIMIT, ANYWINDUP, ANYFILT>
{
    /// Set the PID controller to handle wrapping around between two values, such as `-PI` and `PI`.
    pub fn set_wrapping(self, min: F, max: F) -> Pid<F, Wrapping<F>, ANYLIMIT, ANYWINDUP, ANYFILT> {
        assert!(min < max, "The `min` value must be less than `max`");
        Pid { wrapping: Wrapping { min, max }, ..self }
    }

    /// Set the PID controller to limit its output. Use this in conjunction with `.set_anti_windup(...)`.
    pub fn set_output_limit(self, min: F, max: F) -> Pid<F, ANYWRAP, OutputLimit<F>, ANYWINDUP, ANYFILT> {
        assert!(min < max, "The `min` value must be less than `max`");
        Pid { output_limit: OutputLimit { min, max }, ..self }
    }

    /// Set the PID controller to use an `AntiWindup` method for the integrator.
    pub fn set_anti_windup(self, anti_windup: AntiWindup<F>) -> Pid<F, ANYWRAP, ANYLIMIT, AntiWindup<F>, ANYFILT> {
        match anti_windup {
            AntiWindup::Conditional((min,max)) => assert!(min < max, "The `min` value must be less than `max`")
        }
        Pid { anti_windup, ..self }
    }

    /// Set the PID controller to use a low-pass filter for the derivative term.
    pub fn set_lp_filter(self, tau: F) -> Pid<F, ANYWRAP, ANYLIMIT, ANYWINDUP, LpFilter<F>> {
        assert!(tau > F::zero(), "The time constant `tau` must be strictly positive");
        Pid { lp_filter: LpFilter { tau, signal: Default::default() }, ..self }
    }

    /// Change the gains of the individual controllers.
    pub fn set_gains(&mut self, kp: F, ki: F, kd: F) {
        self.kp = kp;
        self.ki = ki;
        self.kd = kd;
    }

    /// Reset the integral to `0.0`.
    pub fn reset_integral(&mut self) {
        self.integral = Default::default();
    }

    /// Reset the integral to a specific value.
    pub fn reset_integral_to(&mut self, integral: F) {
        self.integral = integral;
    }

    /// Enables/disables the integral.
    pub fn enable_integral(&mut self, enable: bool) {
        self.integral_en = enable;
    }

    /// Enables/disables integrator, and resets integral to `0.0` on disable.
    pub fn enable_reset_integral(&mut self, enable: bool) {
        self.enable_integral(enable);
        self.reset_integral();
    }
}

impl<F: Float, ANYWRAP, ANYLIMIT, ANYWINDUP, ANYFILT> Pid<F, ANYWRAP, ANYLIMIT, ANYWINDUP, ANYFILT>
where
    Pid<F, ANYWRAP, ANYLIMIT, ANYWINDUP, ANYFILT>:
        WrappingTrait<F> + OutputLimitTrait<F> + AntiWindupTrait<F> + LpFilterTrait<F>,
{
    /// Update the controller with an `error` value using a specified sampling time `ts`.
    pub fn update_ts(&mut self, mut error: F, ts: F) -> F {
        // Constrain input
        self.apply_wrapping(&mut error);

        // Proportional gain
        let proportional = self.kp * error;

        // Derivative gain
        let mut derivative = self.kd * (error - self.prev_error) / ts;
        self.apply_lp_filter(&mut derivative, ts);
        self.prev_error = error;

        // Integral gain
        if self.integral_en {
            self.apply_anti_windup(self.ki * error * ts, proportional + derivative);
        }

        // Constrain output
        self.apply_output_limit(proportional + self.integral + derivative)
    }

    /// Update the controller with an `error` value.
    pub fn update(&mut self, error: F) -> F {
        self.update_ts(error, self.ts)
    }
}

// Wrapping of input
pub trait WrappingTrait<F> {
    fn apply_wrapping(&mut self, input: &mut F);
}

impl<F: Float + AddAssign + SubAssign, ANYLIMIT, ANYWINDUP, ANYFILT> WrappingTrait<F>
    for Pid<F, Wrapping<F>, ANYLIMIT, ANYWINDUP, ANYFILT>
{
    #[inline(always)]
    fn apply_wrapping(&mut self, error: &mut F) {
        let minmaxdiff = self.wrapping.max - self.wrapping.min;

        // Wrap the error itself (aggressively)
        let frac = ((*error-self.wrapping.min)/(self.wrapping.max-self.wrapping.min)).floor();
        *error -= frac*(self.wrapping.max-self.wrapping.min);

        // Wrap previous error and integral
        if *error - self.prev_error < self.wrapping.min {
            self.prev_error -= minmaxdiff;
            self.integral -= minmaxdiff;
        } else if *error - self.prev_error > self.wrapping.max {
            self.prev_error += minmaxdiff;
            self.integral += minmaxdiff;
        }
    }
}

impl<F: Float, ANYLIMIT, ANYWINDUP, ANYFILT> WrappingTrait<F>
    for Pid<F, NoWrapping, ANYLIMIT, ANYWINDUP, ANYFILT>
{
    #[inline(always)]
    fn apply_wrapping(&mut self, _input: &mut F) {}
}

// Constraining of output
pub trait OutputLimitTrait<F> {
    fn apply_output_limit(&mut self, input: F) -> F;
}

impl<F: Float + Ord, ANYWRAP, ANYWINDUP, ANYFILT> OutputLimitTrait<F>
    for Pid<F, ANYWRAP, OutputLimit<F>, ANYWINDUP, ANYFILT>
{
    #[inline(always)]
    fn apply_output_limit(&mut self, input: F) -> F {
        input.clamp(self.output_limit.min, self.output_limit.max)
    }
}

impl<F: Float, ANYWRAP, ANYWINDUP, ANYFILT> OutputLimitTrait<F>
    for Pid<F, ANYWRAP, NoOutputLimit, ANYWINDUP, ANYFILT>
{
    #[inline(always)]
    fn apply_output_limit(&mut self, input: F) -> F {
        input
    }
}

// Anti-windup for integral
pub trait AntiWindupTrait<F> {
    fn apply_anti_windup(&mut self, input: F, pd: F);
}

impl<F: Float + AddAssign, ANYWRAP, ANYLIMIT, ANYFILT> AntiWindupTrait<F>
    for Pid<F, ANYWRAP, ANYLIMIT, AntiWindup<F>, ANYFILT>
{
    #[inline(always)]
    fn apply_anti_windup(&mut self, input: F, pd: F) {
        match self.anti_windup {
            AntiWindup::Conditional((min, max)) => {
                if pd > min && pd < max {
                    self.integral += input;
                }
            }
        }
    }
}

impl<F: Float, ANYWRAP, ANYLIMIT, ANYFILT> AntiWindupTrait<F>
    for Pid<F, ANYWRAP, ANYLIMIT, NoAntiWindup, ANYFILT>
{
    #[inline(always)]
    fn apply_anti_windup(&mut self, input: F, _pd: F) {
        self.integral = self.integral + input;
    }
}

// Low pass filter for derivative
pub trait LpFilterTrait<F> {
    fn apply_lp_filter(&mut self, input: &mut F, ts: F);
}

impl<F: Float, ANYWRAP, ANYLIMIT, ANYWINDUP> LpFilterTrait<F>
    for Pid<F, ANYWRAP, ANYLIMIT, ANYWINDUP, LpFilter<F>>
{
    #[inline(always)]
    fn apply_lp_filter(&mut self, input: &mut F, ts: F) {
        let a = ts / (ts + self.lp_filter.tau);
        self.lp_filter.signal = self.lp_filter.signal * (F::one() - a) + *input * a;
        *input = self.lp_filter.signal;
    }
}

impl<F: Float, ANYWRAP, ANYLIMIT, ANYWINDUP> LpFilterTrait<F>
    for Pid<F, ANYWRAP, ANYLIMIT, ANYWINDUP, NoLpFilter>
{
    #[inline(always)]
    fn apply_lp_filter(&mut self, _input: &mut F, _ts: F) {}
}
