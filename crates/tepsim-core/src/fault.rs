//! The twenty disturbances, as what they actually are.
//!
//! `IDV(1..20)` is a bare integer array in the original, and the header at
//! `teprob.f:172-191` names each one in prose. Nothing connects a name to the
//! line that implements it, and five of the twenty are called only "Unknown".
//!
//! This module makes the connection explicit: every fault names the
//! `teprob.f` line it acts on, where it is injected, and what shape it has.
//!
//! # Three shapes, and the third is not a plant disturbance at all
//!
//! [`Shape::Step`] faults change a feed condition the moment they are switched
//! on and hold it: seven of the twenty, at `teprob.f:407-414`, `567` and
//! `568`.
//!
//! [`Shape::Random`] faults enable a walk channel through `IDVWLK`
//! (`teprob.f:347-358`). Ten of the twenty, of which the last three drive
//! *spike trains* rather than walks; see [`mod@crate::walk`].
//!
//! [`Shape::Sticking`] faults do not touch the model. They set `IVST`
//! (`teprob.f:793-798`), which widens the dead band a valve command must cross
//! before the valve follows it. Three of the twenty.
//!
//! That third kind matters more than its size. A sticking fault is not a
//! disturbance to the *plant*; it is a disturbance to the *controller's
//! authority over* the plant. In an open-loop run, where the command never
//! moves, it does nothing whatever, and a scenario engine that treated it as a
//! plant fault would report an injected disturbance with no effect and look
//! broken.
//!
//! # The five "Unknown" faults are not unknown
//!
//! The header calls `IDV(16)` through `IDV(20)` unknown, and every paper on
//! TEP repeats it. The *source* is perfectly explicit about what they do; only
//! their physical interpretation was withheld:
//!
//! | Fault | What `teprob.f` does with it |
//! |---|---|
//! | `IDV(16)` | walk channel 9, the stripper steam valve capacity |
//! | `IDV(17)` | spike channel 10, the reactor coolant duty |
//! | `IDV(18)` | spike channel 11, the condenser coolant duty |
//! | `IDV(19)` | sticks valves 5, 7, 8 and 9 |
//! | `IDV(20)` | spike channel 12, the reactor outlet flow |
//!
//! So `IDV(19)` is a sticking fault and the other four are not, which the
//! shared label hides. Three of the four are the *spike* channels, which is
//! why they are reported in the literature as the hardest to detect: they are
//! intermittent rather than sustained.

extern crate alloc;

/// How a fault enters the model.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Shape {
    /// Changes a feed condition and holds it. `teprob.f:407-414`, `567`,
    /// `568`.
    Step,
    /// Enables a walk channel, which then wanders on its own schedule.
    /// `teprob.f:347-358`.
    Random {
        /// Which channels, one-based. `IDV(8)` and `IDV(13)` drive two each.
        channels: &'static [usize],
        /// Whether those channels are spike trains rather than walks.
        spiking: bool,
    },
    /// Widens a valve's dead band. Touches no equation in the model.
    /// `teprob.f:793-798`.
    Sticking {
        /// Which valves, one-based.
        valves: &'static [usize],
    },
}

/// One of the twenty disturbances.
#[derive(Clone, Copy, Debug)]
pub struct Fault {
    /// The `IDV` index, one-based.
    pub index: usize,
    /// The description from `teprob.f:172-191`, verbatim.
    pub published: &'static str,
    /// What it does, where the header says "Unknown" or where the prose is
    /// less specific than the source.
    pub effect: &'static str,
    /// How it enters.
    pub shape: Shape,
    /// The `teprob.f` line it acts on.
    pub line: &'static str,
}

impl Fault {
    /// Whether this fault reaches the plant model at all.
    ///
    /// False for the three sticking faults: in an open-loop run they do
    /// nothing, and that is not a bug in the scenario.
    #[must_use]
    pub const fn affects_the_plant(&self) -> bool {
        !matches!(self.shape, Shape::Sticking { .. })
    }
}

/// All twenty, in `IDV` order.
///
/// The `published` column is `teprob.f:172-191` verbatim, including the five
/// that say "Unknown"; `effect` is what the source actually does.
//
// @port teprob.f:170-191
pub const FAULTS: [Fault; 20] = [
    Fault {
        index: 1,
        published: "A/C Feed Ratio, B Composition Constant (Stream 4)",
        effect: "steps the mixed feed's A fraction down by 0.03",
        shape: Shape::Step,
        line: "teprob.f:407",
    },
    Fault {
        index: 2,
        published: "B Composition, A/C Ratio Constant (Stream 4)",
        effect: "steps B up by 0.005 and A down by 2.43719e-3, on two lines",
        shape: Shape::Step,
        line: "teprob.f:408-409",
    },
    Fault {
        index: 3,
        published: "D Feed Temperature (Stream 2)",
        effect: "steps the D feed temperature up by 5 C",
        shape: Shape::Step,
        line: "teprob.f:411",
    },
    Fault {
        index: 4,
        published: "Reactor Cooling Water Inlet Temperature",
        effect: "steps the reactor coolant inlet up by 5 C",
        shape: Shape::Step,
        line: "teprob.f:413",
    },
    Fault {
        index: 5,
        published: "Condenser Cooling Water Inlet Temperature",
        effect: "steps the condenser coolant inlet up by 5 C",
        shape: Shape::Step,
        line: "teprob.f:414",
    },
    Fault {
        index: 6,
        published: "A Feed Loss (Stream 1)",
        effect: "shuts the A feed off entirely, not partially",
        shape: Shape::Step,
        line: "teprob.f:567",
    },
    Fault {
        index: 7,
        published: "C Header Pressure Loss - Reduced Availability (Stream 4)",
        effect: "reduces the mixed feed's capacity by 20%",
        shape: Shape::Step,
        line: "teprob.f:568",
    },
    Fault {
        index: 8,
        published: "A, B, C Feed Composition (Stream 4)",
        effect: "enables two walk channels, on A and on B",
        shape: Shape::Random {
            channels: &[1, 2],
            spiking: false,
        },
        line: "teprob.f:347-348",
    },
    Fault {
        index: 9,
        published: "D Feed Temperature (Stream 2)",
        effect: "enables the D feed temperature walk",
        shape: Shape::Random {
            channels: &[3],
            spiking: false,
        },
        line: "teprob.f:349",
    },
    Fault {
        index: 10,
        published: "C Feed Temperature (Stream 4)",
        effect: "enables the mixed feed temperature walk",
        shape: Shape::Random {
            channels: &[4],
            spiking: false,
        },
        line: "teprob.f:350",
    },
    Fault {
        index: 11,
        published: "Reactor Cooling Water Inlet Temperature",
        effect: "enables the reactor coolant inlet walk",
        shape: Shape::Random {
            channels: &[5],
            spiking: false,
        },
        line: "teprob.f:351",
    },
    Fault {
        index: 12,
        published: "Condenser Cooling Water Inlet Temperature",
        effect: "enables the condenser coolant inlet walk",
        shape: Shape::Random {
            channels: &[6],
            spiking: false,
        },
        line: "teprob.f:352",
    },
    Fault {
        index: 13,
        published: "Reaction Kinetics",
        effect: "enables two walks, one per rate constant of reactions 1 and 2",
        shape: Shape::Random {
            channels: &[7, 8],
            spiking: false,
        },
        line: "teprob.f:353-354",
    },
    Fault {
        index: 14,
        published: "Reactor Cooling Water Valve",
        effect: "sticks valve 10; touches no equation in the model",
        shape: Shape::Sticking { valves: &[10] },
        line: "teprob.f:793",
    },
    Fault {
        index: 15,
        published: "Condenser Cooling Water Valve",
        effect: "sticks valve 11; touches no equation in the model",
        shape: Shape::Sticking { valves: &[11] },
        line: "teprob.f:794",
    },
    Fault {
        index: 16,
        published: "Unknown",
        effect: "enables walk channel 9, the stripper steam valve capacity",
        shape: Shape::Random {
            channels: &[9],
            spiking: false,
        },
        line: "teprob.f:355",
    },
    Fault {
        index: 17,
        published: "Unknown",
        effect: "enables spike channel 10, the reactor coolant duty",
        shape: Shape::Random {
            channels: &[10],
            spiking: true,
        },
        line: "teprob.f:356",
    },
    Fault {
        index: 18,
        published: "Unknown",
        effect: "enables spike channel 11, the condenser coolant duty",
        shape: Shape::Random {
            channels: &[11],
            spiking: true,
        },
        line: "teprob.f:357",
    },
    Fault {
        index: 19,
        published: "Unknown",
        effect: "sticks valves 5, 7, 8 and 9; touches no equation in the model",
        shape: Shape::Sticking {
            valves: &[5, 7, 8, 9],
        },
        line: "teprob.f:795-798",
    },
    Fault {
        index: 20,
        published: "Unknown",
        effect: "enables spike channel 12, the reactor outlet flow",
        shape: Shape::Random {
            channels: &[12],
            spiking: true,
        },
        line: "teprob.f:358",
    },
];

/// Look one up by its `IDV` index, one-based.
#[must_use]
pub const fn fault(index: usize) -> Option<&'static Fault> {
    if index >= 1 && index <= FAULTS.len() {
        Some(&FAULTS[index - 1])
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::walk::{CHANNELS, channel_flags};

    /// The table's channel mapping must agree with the code that implements
    /// it. Two statements of one fact, which is the point.
    #[test]
    fn the_table_agrees_with_the_channel_mapping() {
        for entry in &FAULTS {
            let mut idv = [0.0; 20];
            idv[entry.index - 1] = 1.0;
            let flags = channel_flags(&idv);

            let expected: &[usize] = match entry.shape {
                Shape::Random { channels, .. } => channels,
                _ => &[],
            };
            for channel in 1..=CHANNELS {
                let on = flags[channel - 1] == 1;
                assert_eq!(
                    on,
                    expected.contains(&channel),
                    "IDV({}) and channel {channel}: the table says {}, \
                     channel_flags says {on}",
                    entry.index,
                    expected.contains(&channel)
                );
            }
        }
    }

    /// Every channel is driven by exactly one fault, and every fault that
    /// claims a channel gets one. No channel is orphaned.
    #[test]
    fn the_twelve_channels_are_covered_exactly_once() {
        let mut owners = [0_usize; CHANNELS];
        for entry in &FAULTS {
            if let Shape::Random { channels, .. } = entry.shape {
                for channel in channels {
                    owners[channel - 1] += 1;
                }
            }
        }
        for (index, count) in owners.iter().enumerate() {
            assert_eq!(
                *count,
                1,
                "channel {} is driven by {count} faults",
                index + 1
            );
        }
    }

    /// The three spike channels are 10, 11 and 12, and the table says so.
    #[test]
    fn only_the_last_three_channels_spike() {
        for entry in &FAULTS {
            if let Shape::Random { channels, spiking } = entry.shape {
                let all_high = channels.iter().all(|c| *c >= 10);
                assert_eq!(
                    spiking, all_high,
                    "IDV({}) claims spiking = {spiking} for channels {channels:?}",
                    entry.index
                );
            }
        }
    }

    /// The sticking faults touch no equation, and the others all do.
    #[test]
    fn exactly_three_faults_do_not_reach_the_plant() {
        let inert: alloc::vec::Vec<usize> = FAULTS
            .iter()
            .filter(|f| !f.affects_the_plant())
            .map(|f| f.index)
            .collect();
        assert_eq!(inert, alloc::vec![14, 15, 19]);
    }

    /// The five the header calls "Unknown" all have a stated effect here, and
    /// they are not all the same kind.
    #[test]
    fn the_unknown_faults_are_explained_and_are_not_alike() {
        let unknown: alloc::vec::Vec<&Fault> =
            FAULTS.iter().filter(|f| f.published == "Unknown").collect();
        assert_eq!(unknown.len(), 5);
        for entry in &unknown {
            assert!(
                !entry.effect.is_empty() && entry.effect != "Unknown",
                "IDV({}) is still unexplained",
                entry.index
            );
        }
        // One of the five is a sticking fault and the rest are not, which is
        // the distinction the shared label hides.
        let sticking = unknown.iter().filter(|f| !f.affects_the_plant()).count();
        assert_eq!(sticking, 1, "IDV(19) is the odd one out");
        // And three of the five drive spike channels.
        let spiking = unknown
            .iter()
            .filter(|f| matches!(f.shape, Shape::Random { spiking: true, .. }))
            .count();
        assert_eq!(spiking, 3);
    }

    /// Every entry cites a source line, or the table decays into prose nobody
    /// can check.
    #[test]
    fn every_fault_names_the_line_it_acts_on() {
        for (offset, entry) in FAULTS.iter().enumerate() {
            assert_eq!(entry.index, offset + 1, "the table is out of order");
            assert!(
                entry.line.starts_with("teprob.f:"),
                "IDV({}) does not name its line",
                entry.index
            );
            assert!(fault(entry.index).is_some());
        }
        assert!(fault(0).is_none() && fault(21).is_none());
    }
}
