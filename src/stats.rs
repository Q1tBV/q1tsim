extern crate statrs;

use ::std::f64::consts::SQRT_2;

/// Binomial quantile
///
/// Compute lower and upper bounds, such that the chance of finding a value
/// outside these bounds in a binomial distribution of size `nr_shots` and
/// success probability `p`, is less than `tol`.
fn get_bounds(nr_shots: usize, p: f64, tol: f64) -> (usize, usize)
{
    let mu = nr_shots as f64 * p;
    let sigma = (nr_shots as f64 * p * (1.0-p)).sqrt();
    let quantile = mu
        + sigma * SQRT_2 * statrs::function::erf::erf_inv(2.0*tol - 1.0);
    let low = quantile.floor() as usize;
    let high = (mu + (mu - quantile)).ceil() as usize;
    (low, high)
}

/// Check a measurement result
///
/// Check if the measurement result `count` is probably correct, with an error
/// chance of less than `tol`, assuming a binomial distribution of size
/// `nr_shots` and success probability `p`.
/// NOTE: this is a statistical result. It can be expected that tests using this
/// method will fail approximately with a rate of approcimately `tol`.
pub fn measurement_ok(count: usize, nr_shots: usize, p: f64, tol: f64) -> bool
{
    let (low, high) = get_bounds(nr_shots, p, tol);
    count > low && count < high
}

#[cfg(test)]
mod tests
{
    use super::{get_bounds, measurement_ok};

    #[test]
    fn test_get_bounds()
    {
        assert_eq!(get_bounds(1024, 0.5, 1.0e-5), (443, 581));
        assert_eq!(get_bounds(1234, 0.25, 1.0e-5), (243, 374));
        assert_eq!(get_bounds(1234, 0.75, 1.0e-5), (860, 991));
        assert_eq!(get_bounds(1234, 0.75, 1.0e-10), (828, 1023));
        assert_eq!(get_bounds(1_000_000, 0.43, 1.0e-4), (428158, 431842));
    }

    #[test]
    fn test_measurement_ok()
    {
        assert!(measurement_ok(512, 1024, 0.5, 1.0e-5));
        assert!(measurement_ok(444, 1024, 0.5, 1.0e-5));
        assert!(!measurement_ok(443, 1024, 0.5, 1.0e-5));
        assert!(measurement_ok(580, 1024, 0.5, 1.0e-5));
        assert!(!measurement_ok(581, 1024, 0.5, 1.0e-5));
        assert!(!measurement_ok(1023, 1024, 0.5, 1.0e-5));
        assert!(!measurement_ok(0, 1024, 0.5, 1.0e-5));
    }
}
