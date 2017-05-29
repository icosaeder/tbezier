# tbezier
Smooth interpolation without misplaced extremes based on the finite discrete point set

It is an implementation of finite discrete point set smooth interpolation algorithms 
based on cubic Bezier curves with control points calculated according to the tangents to the 
angles of polygonal line that is built by the linear interpolation of the input points set.

Two functions are provided: tbezierSO1 that builds the curve with smoothness order 1 and
tbezierSO0 that builds curve with smoothness order 0 and uses special heuristics to
reduce lengths of tangents and therefore reduce the difference with linear interpolation
making result curve look nicer.

tbezierSO1 is recommended for scientific visualization as it uses strict math to balance
between smoothness and interpolation accuracy.
tbezierSO0 is recommended for advertising purposes as it produces nicer looking curves while
the accuracy is in common case lower.

Read this for algorithm details: http://sv-journal.org/2017-1/04.php?lang=en

Written by Konstantin Ryabinin under terms of MIT license.
