### Version 0.2.0

#### New features

* [feat](https://github.com/EC-USGS/pywatershed/commit/bf763b7ab21f9f0641155e111e1c081d5606f3c3): Automatic release to PyPi workflow. https://github.com/EC-USGS/pywatershed/pull/179. Committed by James McCreight on 2023-05-25.

#### Bug fixes

* [fix](https://github.com/EC-USGS/pywatershed/commit/9c92dabb9aeca2c0281bd9c6e3ac2742ff3e1526): Protect parameters as readonly np.arrays and MappingProxyTypes instead of dicts. Committed by James McCreight on 2023-05-26.

#### Refactoring

* [refactor](https://github.com/EC-USGS/pywatershed/commit/c82b98a528b33b6e81cf29a8234b6bf13c612e85): Refactor dependencies for standard installation. Committed by Joseph Hughes on 2023-05-08.
* [refactor](https://github.com/EC-USGS/pywatershed/commit/e01099277e8e42bd8a5b900eca3ad9170debf910): PRMSChannel init self._muskingum_mann from numpy, numba, fortran instead if during calculate. Committed by James McCreight on 2023-05-26.

#### Self-driving

* [self-driving](https://github.com/EC-USGS/pywatershed/commit/9848b905c85bdba11fc41f54f3e93ab3e65da01e): All but channel budget. Committed by James McCreight on 2023-05-09.
* [self-driving](https://github.com/EC-USGS/pywatershed/commit/33cdd89375640960ad7cd69796ddda16061b40fb): Channel variables to track inputs to channel and fix mass balance. Committed by James McCreight on 2023-05-09.

