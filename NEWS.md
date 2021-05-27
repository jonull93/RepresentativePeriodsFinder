RepresentativePeriodsFinder release notes
==========================================

Version 0.4.2 (May 27th 2021)
- Fixed ordering error bug
- Added example for re-ordering

Version 0.4.1 (March 24th 2021)
---------------------------------
- There was an error in the optimisation formulation when re-ordering, this was fixed.
- Fixed deployment of documentation.
- Restricted versions of DataFrames and CSV since these can cause issues. Also added entries about this in README.md
- Extended tests to be performed also for Julia 1.2.0, 1.3.0 and 1.4.0 (but only on the master branch).
- Changed `Getting Started` in documentation to allow for copy pasting the `.yaml` configuration file.

Version 0.4.0 (February 15th 2021)
---------------------------------
- General overhaul of code and increased documentation.