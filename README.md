# gfsieve
An OpenCL implementation of a Generalized Fermat Number sieving program

## About

**gfsieve** is an [OpenCL™](https://www.khronos.org/opencl/) application.  

It finds the prime factors of the form *b*<sup>2<sup>n</sup></sup>&nbsp;+&nbsp;1, for even 2&nbsp;&le;&nbsp;*b*&nbsp;&le;&nbsp;2&thinsp;&middot;&thinsp;10<sup>9</sup> and fixed 15&nbsp;&le;&nbsp;*n*&nbsp;&le;&nbsp;24.  
Prime factors *p* are in the range 10<sup>15</sup>&nbsp;&le;&nbsp;*p*&nbsp;&le;&nbsp;6&thinsp;&middot;&thinsp;10<sup>23</sup>.  
A CPU application is available for *p*&nbsp;&le;&nbsp;10<sup>15</sup>.  

It is based on Phil Carmody's algorithm (http://fatphil.org/maths/GFN/maths.html).  

gfsieve is a highly optimised GPU application, created in 2020.

## Build

This version was compiled with gcc and tested on Windows and Linux (Ubuntu).  

An OpenCL SDK is not required. OpenCL header files are included in the project and the application can be linked with the dynamic OpenCL library of the OS.

## TODO
