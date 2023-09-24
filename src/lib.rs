#![feature(type_changing_struct_update)]
#![no_std]

pub struct NoCircular;
pub struct Circular {
    pub min: f32,
    pub max: f32,
}

pub struct NoOutputLimit;
pub struct OutputLimit {
    pub min: f32,
    pub max: f32,
}

pub struct NoAntiWindup;
pub enum AntiWindup {
    Conditional((f32, f32)),
}

pub struct NoLpFilter;
pub struct LpFilter {
    tau : f32,
    signal : f32,
}

pub struct Pid<CIRC,LIMIT,WINDUP,LPFILT> {
    kp: f32,
    ki: f32,
    kd: f32,
    ts: f32,
    integral: f32,
    integral_en: bool,
    prev_error: f32,
    circular: CIRC,
    output_limit: LIMIT,
    anti_windup: WINDUP,
    lp_filter : LPFILT,
}

impl Pid<NoCircular,NoOutputLimit,NoAntiWindup,NoLpFilter> {
    pub fn new(kp:f32, ki:f32, kd:f32, ideal:bool, ts:f32) -> Pid<NoCircular,NoOutputLimit,NoAntiWindup,NoLpFilter> {
        Pid {
            kp,
            ki : if ideal {ki*kp} else {ki},
            kd : if ideal {kd*kp} else {kd} ,
            ts,
            integral : 0.,
            integral_en : true,
            prev_error : 0.,
            circular : NoCircular,
            output_limit : NoOutputLimit,
            anti_windup : NoAntiWindup,
            lp_filter : NoLpFilter
        }
    }
}

impl<ANYCIRC,ANYLIMIT,ANYWINDUP,ANYFILT> Pid<ANYCIRC,ANYLIMIT,ANYWINDUP,ANYFILT> {


    pub fn set_circular(self, min: f32, max : f32) -> Pid<Circular,ANYLIMIT,ANYWINDUP,ANYFILT> {
        Pid{ circular : Circular { min, max } , ..self }
    }

    pub fn set_output_limit(self, output_limit: OutputLimit) -> Pid<ANYCIRC,OutputLimit,ANYWINDUP,ANYFILT> {
        Pid{ output_limit , ..self }
    }

    pub fn set_anti_windup(self, anti_windup: AntiWindup) -> Pid<ANYCIRC,ANYLIMIT,AntiWindup,ANYFILT> {
        Pid{ anti_windup , ..self }
    }

    pub fn set_lp_filter(self, tau: f32) -> Pid<ANYCIRC,ANYLIMIT,ANYWINDUP,LpFilter> {
        Pid{ lp_filter : LpFilter{ tau, signal : 0. } , ..self }
    }

    pub fn set_gains(&mut self, kp: f32, ki: f32, kd: f32) {
        self.kp = kp;
        self.ki = ki;
        self.kd = kd;
    }

    pub fn reset_integral(& mut self) {
        self.integral = 0.0;
    }

    pub fn reset_integral_to(& mut self, integral : f32) {
        self.integral = integral;
    }

    pub fn enable_integral(& mut self, enable : bool) {
        self.integral_en = enable;
    }
}

impl <ANYCIRC,ANYLIMIT,ANYWINDUP,ANYFILT> Pid<ANYCIRC,ANYLIMIT,ANYWINDUP,ANYFILT>
where
    Pid<ANYCIRC,ANYLIMIT,ANYWINDUP,ANYFILT> : CircularTrait + OutputLimitTrait + AntiWindupTrait + LpFilterTrait
{

    pub fn update_ts(&mut self, error: f32, ts: f32) -> f32 {

        // Constrain input
        self.apply_circular(&error);

        // Proportional gain
        let proportional = self.kp * error;

        // Derivative gain
        let mut derivative =  self.kd * (error - self.prev_error) / ts;
        self.apply_lp_filter(& mut derivative, &ts);
        self.prev_error = error;

        // Integral gain
        if self.integral_en {
            self.apply_anti_windup( self.ki * error * ts , proportional + derivative);
        }

        // Constrain output
        self.apply_output_limit(proportional + self.integral + derivative)

    }

    pub fn update(&mut self, error: f32) -> f32 {

       self.update_ts(error, self.ts)

    }
}

// Wrapping of input
pub trait CircularTrait {
    fn apply_circular(&mut self, input: &f32);
}

impl<ANYLIMIT,ANYWINDUP,ANYFILT> CircularTrait for Pid<Circular,ANYLIMIT,ANYWINDUP,ANYFILT> {
    #[inline(always)]
    fn apply_circular(&mut self, error: &f32) {
        if error - self.prev_error < self.circular.min {
            self.prev_error -= self.circular.max - self.circular.min;
            self.integral -= self.circular.max - self.circular.min;
        } else if error - self.prev_error > self.circular.max {
            self.prev_error += self.circular.max - self.circular.min;
            self.integral += self.circular.max - self.circular.min;
        }
    }
}

impl<ANYLIMIT,ANYWINDUP,ANYFILT> CircularTrait for Pid<NoCircular,ANYLIMIT,ANYWINDUP,ANYFILT> {

    #[inline(always)]
    fn apply_circular(&mut self, _input: &f32) {}
}

// Constraining of output
pub trait OutputLimitTrait {
    fn apply_output_limit(&mut self, input: f32 ) -> f32;
}

impl<ANYCIRC,ANYWINDUP,ANYFILT> OutputLimitTrait for Pid<ANYCIRC,OutputLimit,ANYWINDUP,ANYFILT> {
    #[inline(always)]
    fn apply_output_limit(&mut self, input: f32 ) -> f32 {
        input.clamp(self.output_limit.min, self.output_limit.max)
    }
}

impl<ANYCIRC,ANYWINDUP,ANYFILT> OutputLimitTrait for Pid<ANYCIRC,NoOutputLimit,ANYWINDUP,ANYFILT> {

    #[inline(always)]
    fn apply_output_limit(&mut self, input: f32  ) -> f32 {
        input
    }
}

// Anti-windup for integral
pub trait AntiWindupTrait {
    fn apply_anti_windup(&mut self, input: f32, pd : f32 );
}

impl<ANYCIRC,ANYLIMIT,ANYFILT> AntiWindupTrait for Pid<ANYCIRC,ANYLIMIT,AntiWindup,ANYFILT> {
    #[inline(always)]
    fn apply_anti_windup(&mut self, input: f32 , pd : f32 ) {
        match self.anti_windup {
            AntiWindup::Conditional((min,max)) => {
                if pd > min && pd < max {
                    self.integral += input;
                }
            },
        }
    }
}

impl<ANYCIRC,ANYLIMIT,ANYFILT> AntiWindupTrait for Pid<ANYCIRC,ANYLIMIT,NoAntiWindup,ANYFILT> {

    #[inline(always)]
    fn apply_anti_windup(&mut self, input: f32 , _pd : f32 ) {
        self.integral += input;
    }
}

// Low pass filter for derivative
pub trait LpFilterTrait {
    fn apply_lp_filter(&mut self, input: &mut f32, ts: &f32);
}

impl<ANYCIRC,ANYLIMIT,ANYWINDUP> LpFilterTrait for Pid<ANYCIRC,ANYLIMIT,ANYWINDUP,LpFilter> {
    #[inline(always)]
    fn apply_lp_filter(&mut self, input: &mut f32, ts: &f32) {

        let a = ts/(ts+self.lp_filter.tau);
        self.lp_filter.signal = self.lp_filter.signal*(1.0-a) + *input*a;
        *input = self.lp_filter.signal;
    }
}

impl<ANYCIRC,ANYLIMIT,ANYWINDUP> LpFilterTrait for Pid<ANYCIRC,ANYLIMIT,ANYWINDUP,NoLpFilter> {

    #[inline(always)]
    fn apply_lp_filter(&mut self, _input: &mut f32, _ts: &f32) {}
}
