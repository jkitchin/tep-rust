//! The walk spans and noise magnitudes against `COMMON/WLK/` and
//! `COMMON/TEPROC/` as gfortran actually filled them.
//!
//! `tepsim-core/tests/walk_constants.rs` reparses the Fortran and derives what
//! each value should be. This one asks the compiler. It is the decisive check
//! of the two, because it does not share the other's reading of Fortran
//! literal semantics: if both are wrong in the same way, only this one
//! notices.

#![cfg(feature = "oracle")]

use tepsim_core::CHANNEL_SPANS;
use tepsim_core::constants::MEASUREMENT_NOISE;
use tepsim_oracle::Oracle;

#[test]
fn the_noise_magnitudes_are_what_gfortran_stored() {
    let mut oracle = Oracle::lock();
    // `TEINIT` is what fills them; a bare lock has zeros.
    let _ = oracle.init();
    let theirs = oracle.teproc().xns;

    for index in 0..41 {
        assert_eq!(
            MEASUREMENT_NOISE[index].to_bits(),
            theirs[index].to_bits(),
            "XNS({}): ours {}, gfortran {}",
            index + 1,
            MEASUREMENT_NOISE[index],
            theirs[index]
        );
    }
    assert!(
        theirs.iter().any(|x| *x > 0.0),
        "every magnitude is zero, so TEINIT did not run and this proves nothing"
    );
}

#[test]
fn the_channel_spans_are_what_gfortran_stored() {
    let mut oracle = Oracle::lock();
    let _ = oracle.init();
    let wlk = oracle.wlk();

    for channel in 0..12 {
        let ours = &CHANNEL_SPANS[channel];
        let cases = [
            ("HSPAN", ours.duration_span, wlk.hspan[channel]),
            ("HZERO", ours.duration_centre, wlk.hzero[channel]),
            ("SSPAN", ours.value_span, wlk.sspan[channel]),
            ("SZERO", ours.value_centre, wlk.szero[channel]),
            ("SPSPAN", ours.slope_span, wlk.spspan[channel]),
        ];
        for (name, ours, theirs) in cases {
            assert_eq!(
                ours.to_bits(),
                theirs.to_bits(),
                "{name}({}): ours {ours}, gfortran {theirs}",
                channel + 1
            );
        }
    }
    assert!(
        wlk.hzero.iter().all(|h| *h > 0.0),
        "the spans are still zero, so TEINIT did not run"
    );
}

/// The initial walk state: every channel starts at its own `SZERO`, flat, with
/// its first segment ending at 0.1 hours.
///
/// `teprob.f:1360-1367`. It is not in the constant tables because it is the
/// *state*, and B-0031 will carry it; asserting it here pins what B-0031 has
/// to reproduce.
#[test]
fn the_initial_walk_state_is_flat_at_each_channels_centre() {
    let mut oracle = Oracle::lock();
    let _ = oracle.init();
    let wlk = oracle.wlk();

    for channel in 0..12 {
        assert_eq!(
            wlk.adist[channel].to_bits(),
            CHANNEL_SPANS[channel].value_centre.to_bits(),
            "ADIST({}) does not start at SZERO",
            channel + 1
        );
        for (name, value) in [
            ("BDIST", wlk.bdist[channel]),
            ("CDIST", wlk.cdist[channel]),
            ("DDIST", wlk.ddist[channel]),
            ("TLAST", wlk.tlast[channel]),
        ] {
            assert_eq!(
                value.to_bits(),
                0.0_f64.to_bits(),
                "{name}({}) starts at {value}, not zero",
                channel + 1
            );
        }
        assert_eq!(
            wlk.tnext[channel].to_bits(),
            0.1_f64.to_bits(),
            "TNEXT({}) does not start at 0.1",
            channel + 1
        );
    }
}

/// `IDVWLK` starts at zero for every channel, so a fresh plant has no
/// disturbance active and every walk is a flat line at its centre.
///
/// That is what makes a `d00` run a nominal run.
#[test]
fn a_fresh_plant_has_no_disturbance_active() {
    let mut oracle = Oracle::lock();
    let _ = oracle.init();
    let wlk = oracle.wlk();
    assert!(
        wlk.idvwlk.iter().all(|f| *f == 0),
        "IDVWLK starts as {:?}, so some channel is disturbed before anything \
         asked for it",
        wlk.idvwlk
    );
    assert!(
        oracle.disturbances().iter().all(|d| *d == 0),
        "IDV starts non-zero"
    );
}
