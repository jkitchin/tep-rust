//! The fault table's claims, checked against the Fortran.
//!
//! B-0033. `tepsim_core::fault` says what each of the twenty disturbances
//! does. That is a claim about the original, so it is measured against the
//! original rather than asserted from a reading of it.
//!
//! The header at `teprob.f:172-191` calls five of the twenty "Unknown", and
//! the literature repeats it. The source is explicit about all five; only
//! their physical interpretation was withheld. The tests here demonstrate the
//! difference between the two statements.

#![cfg(feature = "oracle")]

use tepsim_core::fault::{FAULTS, Shape};
use tepsim_core::walk::CHANNELS;
use tepsim_oracle::Oracle;
use tepsim_oracle::tier2::{Pools, Scenario};

const DT: f64 = 1.0 / 3600.0;

/// Force a scenario with one fault on and report what moved.
fn with_fault(
    oracle: &mut Oracle,
    base: &Scenario,
    fault: usize,
) -> (Vec<f64>, [i32; 12], [i32; 12]) {
    let mut scenario = base.clone();
    scenario.disturbances[fault - 1] = 1;
    let snapshot = scenario.force(oracle);
    let wlk = oracle.wlk();
    (
        snapshot.derivative.to_vec(),
        wlk.idvwlk,
        oracle.teproc().ivst,
    )
}

/// Every fault's claimed channels are the channels `IDVWLK` actually enables.
#[test]
fn the_claimed_channels_are_the_ones_the_fortran_enables() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 10, DT);
    let base = pools.nominal_case(0);

    for entry in &FAULTS {
        let (_, idvwlk, _) = with_fault(&mut oracle, &base, entry.index);
        let expected: &[usize] = match entry.shape {
            Shape::Random { channels, .. } => channels,
            _ => &[],
        };
        for channel in 1..=CHANNELS {
            assert_eq!(
                idvwlk[channel - 1] == 1,
                expected.contains(&channel),
                "IDV({}) and channel {channel}: the table claims {}, the \
                 Fortran set IDVWLK to {}",
                entry.index,
                expected.contains(&channel),
                idvwlk[channel - 1]
            );
        }
    }
}

/// Every fault's claimed valves are the valves `IVST` actually sticks.
#[test]
fn the_claimed_valves_are_the_ones_the_fortran_sticks() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 10, DT);
    let base = pools.nominal_case(0);

    for entry in &FAULTS {
        let (_, _, ivst) = with_fault(&mut oracle, &base, entry.index);
        let expected: &[usize] = match entry.shape {
            Shape::Sticking { valves } => valves,
            _ => &[],
        };
        for valve in 1..=12 {
            assert_eq!(
                ivst[valve - 1] == 1,
                expected.contains(&valve),
                "IDV({}) and valve {valve}: the table claims {}, the Fortran \
                 set IVST to {}",
                entry.index,
                expected.contains(&valve),
                ivst[valve - 1]
            );
        }
    }
}

/// The three sticking faults change *nothing* in an open-loop run.
///
/// This is the claim worth checking against the Fortran rather than reasoning
/// about, and it is why they get their own shape in the table. A scenario
/// engine that treated them as plant disturbances would report an injected
/// fault with no effect and look broken.
#[test]
fn a_sticking_fault_changes_nothing_when_the_command_never_moves() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 60, DT);
    let base = pools.nominal_case(30);
    let clean = base.force(&mut oracle).derivative;

    for entry in FAULTS.iter().filter(|f| !f.affects_the_plant()) {
        let (derivative, _, _) = with_fault(&mut oracle, &base, entry.index);
        for (slot, (a, b)) in derivative.iter().zip(clean).enumerate() {
            assert_eq!(
                a.to_bits(),
                b.to_bits(),
                "IDV({}) moved YP({}) in an open-loop run, so it does reach \
                 the plant after all and the table is wrong",
                entry.index,
                slot + 1
            );
        }
    }

    // And every fault that *is* claimed to reach the plant does move
    // something, over a run long enough for a walk to fire. Otherwise the
    // distinction is untested in the other direction.
    let mut inert = Vec::new();
    for entry in FAULTS.iter().filter(|f| f.affects_the_plant()) {
        let mut scenario = base.clone();
        scenario.disturbances[entry.index - 1] = 1;
        let mut moved = false;
        let mut t = 0.0;
        for _ in 0..400 {
            scenario.time = t;
            let with = scenario.force(&mut oracle).derivative;
            scenario.walk = oracle.wlk();
            scenario.rng = oracle.rng();
            scenario.common = oracle.teproc();
            scenario.measurements = oracle.measurements();

            let mut without = scenario.clone();
            without.disturbances = [0; 20];
            let plain = without.force(&mut oracle).derivative;
            if with
                .iter()
                .zip(plain)
                .any(|(a, b)| a.to_bits() != b.to_bits())
            {
                moved = true;
                break;
            }
            t += 0.05;
        }
        if !moved {
            inert.push(entry.index);
        }
    }
    assert!(
        inert.is_empty(),
        "IDV{inert:?} are claimed to reach the plant and moved nothing in 20 \
         simulated hours"
    );
}

/// The step faults act immediately; the random ones need a channel to come
/// due first.
///
/// That difference is the whole distinction between `Shape::Step` and
/// `Shape::Random`, and it is visible in the very first evaluation.
#[test]
fn step_faults_act_at_once_and_random_ones_do_not() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 60, DT);
    // A time before any channel is due: `TNEXT` starts at 0.1 for all twelve.
    let base = Scenario {
        time: 0.05,
        ..pools.nominal_case(0)
    };
    let clean = base.force(&mut oracle).derivative;

    for entry in &FAULTS {
        let (derivative, _, _) = with_fault(&mut oracle, &base, entry.index);
        let moved = derivative
            .iter()
            .zip(clean)
            .any(|(a, b)| a.to_bits() != b.to_bits());
        match entry.shape {
            Shape::Step => assert!(
                moved,
                "IDV({}) is a step fault and changed nothing immediately",
                entry.index
            ),
            Shape::Random { .. } => assert!(
                !moved,
                "IDV({}) is a random-variation fault and changed something \
                 before any channel came due; a walk flag should only take \
                 effect at the next re-segmentation",
                entry.index
            ),
            Shape::Sticking { .. } => assert!(!moved),
        }
    }
}

/// The five "Unknown" faults are three different kinds, which the shared label
/// hides.
#[test]
fn the_five_unknown_faults_are_not_one_kind() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 10, DT);
    let base = pools.nominal_case(0);

    let mut kinds = std::collections::BTreeSet::new();
    for entry in FAULTS.iter().filter(|f| f.published == "Unknown") {
        let (_, idvwlk, ivst) = with_fault(&mut oracle, &base, entry.index);
        let channels: Vec<usize> = (1..=CHANNELS).filter(|c| idvwlk[c - 1] == 1).collect();
        let valves: Vec<usize> = (1..=12).filter(|v| ivst[v - 1] == 1).collect();
        println!(
            "IDV({}) '{}': channels {channels:?}, valves {valves:?} -- {}",
            entry.index, entry.published, entry.effect
        );
        kinds.insert(match entry.shape {
            Shape::Step => "step",
            Shape::Random { spiking: true, .. } => "spike",
            Shape::Random { .. } => "walk",
            Shape::Sticking { .. } => "sticking",
        });
    }
    assert_eq!(
        kinds.len(),
        3,
        "the five 'Unknown' faults came out as {kinds:?}; the table says walk, \
         spike and sticking"
    );
}
