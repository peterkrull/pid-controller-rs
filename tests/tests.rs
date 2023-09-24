#[cfg(test)]
mod tests {

    use assert_approx_eq::assert_approx_eq;
    use pid_controller_rs::*;

    #[test]
    fn test_proportional() {

        let mut controller = Pid::new(2.5, 0.0, 0.0, false,0.25);

        let output = controller.update_ts(3.0, 0.25);
        assert_approx_eq!(output, 7.5, 1e-6);
    }

    #[test]
    fn test_derivative() {

        let mut controller = Pid::new(2.5, 0.0, 0.5, false,0.25);

        let mut output = controller.update_ts(3.0, 0.25);
        assert_approx_eq!(output, 7.5 + 3.0*0.5/0.25, 1e-6);

        output = controller.update_ts(3.0, 0.25);
        assert_approx_eq!(output, 7.5, 1e-6);

        output = controller.update_ts(0.0, 0.25);
        assert_approx_eq!(output, - 3.0*0.5/0.25, 1e-6);

    }

    #[test]
    fn test_integral() {

        let mut controller = Pid::new(2.5, 0.2, 0.0, false,0.25);

        let mut output = controller.update_ts(3.0, 0.25);
        assert_approx_eq!(output, 7.5 + 3.0*0.2*0.25, 1e-6);

        output = controller.update_ts(3.0, 0.25);
        assert_approx_eq!(output, 7.5 + 6.0*0.2*0.25, 1e-6);

        output = controller.update_ts(0.0, 0.25);
        assert_approx_eq!(output, 6.0*0.2*0.25, 1e-6);

    }

    #[test]
    fn test_integral_ideal() {

        let mut controller = Pid::new(2.5, 0.2, 0.0, true,0.25);

        let mut output = controller.update_ts(3.0, 0.25);
        assert_approx_eq!(output, 7.5 + 2.5*3.0*0.2*0.25, 1e-6);

        output = controller.update_ts(3.0, 0.25);
        assert_approx_eq!(output, 7.5 + 2.5*6.0*0.2*0.25, 1e-6);

        output = controller.update_ts(0.0, 0.25);
        assert_approx_eq!(output, 2.5*6.0*0.2*0.25, 1e-6);

    }

    #[test]
    fn test_integral_ideal_0kp() {

        let mut controller = Pid::new(0.0, 0.2, 5.0, true,0.25);

        let mut output = controller.update_ts(3.0, 0.25);
        assert_approx_eq!(output, 0.0, 1e-6);

        output = controller.update_ts(3.0, 0.25);
        assert_approx_eq!(output, 0.0, 1e-6);

        output = controller.update_ts(0.0, 0.25);
        assert_approx_eq!(output, 0.0, 1e-6);

    }
}
