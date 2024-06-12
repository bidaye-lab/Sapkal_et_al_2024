## Leg kinematic analyses 
### Stop neruon activation in intact flies 
- Input: parquet file with all raw data ; format documented [elsewhere](https://github.com/bidaye-lab/kinematics_analysis/blob/main/docs/data_structure.md)
- Segments the data into trials where optogenetic stimulation onset happened when each leg was in the swing or stance phase (Fig 4 a-b)
- Calculates swing duration pre and post optogenetic stimulation (Fig 4d)
- Plots the Femur-Tibia flexion angle probability density over the entire trial duration (fig 4c, extended fig 7a)
- Calculates stopping bout related parameters (Extended Fig 3 c-e)

### Coativation in decapitated flies
- BRK + MDN coactivation (imports parquet file with raw data, refines step cycles predictions based on joint angles, plots kinematic parameters)
- BRK + BDN2 coactivation (imports parquet file with raw data, plots kinematic parameters)

## Scripts for muscle imaging
- muscle imaging scripts are [here](https://github.com/bidaye-lab/Sapkal_et_al_2024/tree/main/Figure4/Muscle-Imaging)
