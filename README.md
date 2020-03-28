# gfsieve
An OpenCL implementation of a Generalized Fermat Number sieving program

## About

**gfsieve** is an [OpenCL™](https://www.khronos.org/opencl/) application.  
It finds prime factors for a range of numbers of the form *b*<sup>2<sup>n</sup></sup> + 1. It is based on Phil Carmody's algorithm (http://fatphil.org/maths/GFN/maths.html).  
gfsieve is a highly optimised GPU application, created in 2020.

## Build

This version was compiled with gcc and tested on Windows and Linux (Ubuntu).  
An OpenCL SDK is not required. OpenCL header files are included in the project and the application can be linked with the dynamic OpenCL library of the OS.

## TODO
