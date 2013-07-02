#!/bin/tcsh -f
foreach x (*.eps)
 convert $x `basename $x eps`png
end
