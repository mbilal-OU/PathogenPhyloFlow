# CI results preview

This directory contains a compact, committed preview from the successful four-genome continuous-integration subset of the larger 10-genome O157:H7 tutorial. The subset includes EDL933, Sakai, EC4115, and TW14359.

These files are included so the figures shown in the main README are traceable to actual workflow outputs rather than mock results.

Key validation results:

- automatic reference selection: EDL933
- accessory gene families: 5,756
- recombination-masked positions: 144,481 of 5,620,522 alignment positions (2.57%)
- raw versus filtered topology: normalized Robinson-Foulds distance 0
- filtered/raw total branch-length ratio: 0.528
- SNP/accessory pairs evaluated: 6
- pairs crossing both configured discordance thresholds: 0
- temporal diagnostic: R² 0.0016; configured screen not passed

The panel is for workflow validation. It is not a reconstruction of one epidemiologically linked outbreak, and these outputs should not be interpreted as direct transmission evidence.
