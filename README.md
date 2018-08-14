## LAPACK Helper module

This module contains subroutines to perform matrix operations using the LAPACK library.
Distributed under the GNU GENERAL PUBLIC LICENSE.

Author: B. Font Garcia
September 2016

Some procedures within the module have been extracted from other sources.
See the Wiki for more details: https://github.com/b-fg/LAPACK-Helper/wiki

### LAPACK Library installation
 - LAPACK can be installed with make. Configuration have to be set in the
 make.inc file. A make.inc.example for a Linux machine running GNU compilers
 is given in the main directory. Some specific make.inc are also available in
 the INSTALL directory
 - LAPACK includes also the CMAKE build. You will need to have CMAKE installed
 on your machine (CMAKE is available at http://www.cmake.org/). CMAKE will allow
 an easy installation on a Windows Machine
 - Specific information to run LAPACK under Windows are available at
 http://icl.cs.utk.edu/lapack-for-windows/lapack/

 For further information on LAPACK please read our FAQ at
 http://www.netlib.org/lapack/#_faq
 A User forum is also available to help you with the LAPACK library at
 http://icl.cs.utk.edu/lapack-forum/
