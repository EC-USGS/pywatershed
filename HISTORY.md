### Version 0.2.0

#### New features

* [feat](https://github.com/EC-USGS/pywatershed/commit/bf763b7ab21f9f0641155e111e1c081d5606f3c3): Automatic release to PyPi workflow. https://github.com/EC-USGS/pywatershed/pull/179. Committed by James McCreight on 2023-05-25.
* [feat](https://github.com/EC-USGS/pywatershed/commit/3b085e29db2ca8901a355187991e9d6df8084955): Individual process models run from separated parameters.. Committed by James McCreight on 2023-06-07.
* [feat](https://github.com/EC-USGS/pywatershed/commit/916ff976167ffb99dd1e7a45f2df8033b6233611): Model and Control from yaml files. Committed by James McCreight on 2023-06-21.
* [feat](https://github.com/EC-USGS/pywatershed/commit/998efbaf84320dda0f545bac3e8931ff211c5dee): Add pre-commit hook to strip output from jupyter notebooks. Committed by James McCreight on 2023-06-23.
* [feat](https://github.com/EC-USGS/pywatershed/commit/f53a5a96e934e6121bded7fe4a8bf452cd9e63d0): New notebook 01_multi-process_models. Committed by James McCreight on 2023-06-30.

#### Bug fixes

* [fix](https://github.com/EC-USGS/pywatershed/commit/9c92dabb9aeca2c0281bd9c6e3ac2742ff3e1526): Protect parameters as readonly np.arrays and MappingProxyTypes instead of dicts. Committed by James McCreight on 2023-05-26.

#### Reactor

* [reactor](https://github.com/EC-USGS/pywatershed/commit/c42808eb8f6e234f3999eb4fd0b777e77d008b2a): StorageUnit renamed to Process, before going to ConservativeProcess subclass Process is conservative. Committed by James McCreight on 2023-06-22.

#### Refactoring

* [refactor](https://github.com/EC-USGS/pywatershed/commit/c82b98a528b33b6e81cf29a8234b6bf13c612e85): Refactor dependencies for standard installation. Committed by Joseph Hughes on 2023-05-08.
* [refactor](https://github.com/EC-USGS/pywatershed/commit/e01099277e8e42bd8a5b900eca3ad9170debf910): PRMSChannel init self._muskingum_mann from numpy, numba, fortran instead if during calculate. Committed by James McCreight on 2023-05-26.
* [refactor](https://github.com/EC-USGS/pywatershed/commit/c5a0ce10f3bf487d7cd508d7daac9494d92102ce): Separated PRMS parameters into discretizations for hrus and segments (dis_hru, dis_seg) and onto individual PRMS processes defined in pywatershed. Committed by James McCreight on 2023-06-14.
* [refactor](https://github.com/EC-USGS/pywatershed/commit/6be7331a8ec40835d7f0a4061c4b9cf5285b5e68): Control does not take/know parameters, pass dis and parameters to individual processes, individual process tests passing. Committed by James McCreight on 2023-06-15.
* [refactor](https://github.com/EC-USGS/pywatershed/commit/afbe91c296f3f78e211e7fa75e2c2dc81fd5375d): Model now takes a dictionary of [control, dis, process, order] while maintaining backwards compatability, very minor changes to api (passing list instead of unpacking list for arguments).. Committed by James McCreight on 2023-06-16.
* [refactor](https://github.com/EC-USGS/pywatershed/commit/cc44610c2d1c665336fe33074e4bd8d0908c2c0c): Calc_method clean up on channel and gw. Committed by James McCreight on 2023-06-22.
* [refactor](https://github.com/EC-USGS/pywatershed/commit/ff872fabfd406ddcf96653e7434034402c9e00ba): Calc_method clean up on canopy. Committed by James McCreight on 2023-06-22.
* [refactor](https://github.com/EC-USGS/pywatershed/commit/7310f2c7e8194bcfd5465f985dd411cebdd5b7f1): Calc_method clean up on prms snow. Committed by James McCreight on 2023-06-22.
* [refactor](https://github.com/EC-USGS/pywatershed/commit/cedd1e72f9a464e758e98698a514ea874dc63d15): Calc_method clean up on prms runoff and soilzone. Committed by James McCreight on 2023-06-22.
* [refactor](https://github.com/EC-USGS/pywatershed/commit/a78c75ab8fec90501dcb3a1d30fa63df48a23f63): Rename StorageUnit to ConservativeProcess class, subclass it from a Process class which does not have a budget. remove options from model except for input_dir, all other options go via control, new set_options on Process and ConservativeProcess. Committed by James McCreight on 2023-06-23.

#### Self-driving

* [self-driving](https://github.com/EC-USGS/pywatershed/commit/9848b905c85bdba11fc41f54f3e93ab3e65da01e): All but channel budget. Committed by James McCreight on 2023-05-09.
* [self-driving](https://github.com/EC-USGS/pywatershed/commit/33cdd89375640960ad7cd69796ddda16061b40fb): Channel variables to track inputs to channel and fix mass balance. Committed by James McCreight on 2023-05-09.

