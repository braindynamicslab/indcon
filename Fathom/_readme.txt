FATHOM TOOLBOX FOR MATLAB

The 'Fathom Toolbox for Matlab' is a collection of Matlab functions and scripts I've written for my every day work as a fisheries oceanographer and fish ecologist. I am releasing it to the public under the GNU General Public License (version 2) to encourage the sharing of code and prevent duplication of effort.

If you find this toolbox useful and/or use it to produce a report or publication, please drop me a line. I'd also appreciate bug reports and suggestions for improvements. While I've made every attempt to write code that provides accurate and precise results, the functions in this toolbox are provided as is, with no guarantees.

Suggested citation:
Jones, D. L. 2012. The Fathom Toolbox for Matlab: multivariate ecological and oceanographic data analysis. College of Marine Science, University of South Florida, St. Petersburg, Florida, USA. Available from: http://www.marine.usf.edu/user/djones/

To install the toolbox:
1) unzip the compressed FTM.zip file
2) place the 'Fathom' folder on your computer where your Matlab toolboxes are stored
3) Add the 'Fathom' folder to your Matlab path: File > Set Path... > Add Folder...
Use the Matlab 'help' command to learn what arguments are required by each function. For example:

==========================================================================
>> help f_stnd
  - standardize values of a matrix, column-wise (= z-scores)
 
  Usage: Xstd = f_stnd(x,{y});
 
  x    = matrix to standardize
  y    = optional, when present X will be standardized according to Y
 
  Xstd = standardized matrix
 
  SEE ALSO: f_center, f_ranging, f_transform
==========================================================================

Most functions have detailed documentation and references added as text comments at the beginning of each file. To view this information, load the *.m Matlab function file in the Matlab editor or any text editor. Check out the 'examples' folder located within the main 'Fathom' folder to see demonstrations of how to use many of the functions. Additional information is also available in the 'doc' folder.

David L. Jones, PhD
djones14@mail.usf.edu
