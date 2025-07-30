# Frequently asked questions FAQ

## Why should I use Julia?

**Short answer:** Speed (this code is 3-5 times faster han Matlab and at least 10 times faster than Python), which is very useful for solving optimization problems.
**Long answer:** Read [Why am I using Julia](https://ufechner7.github.io/2022/08/13/why-julia.html).

## What is a Julia environment?
A Julia environment is a folder that contains a `Project.toml` and a `Manifest.toml` file.
The `Project.toml` contains a list of the packages that are needed to run your code. It can also contain a list of the compatible package versions and more. It is important to use separate projects for each of your pieces of software that you develop, if you are not doing that (and use only the global environment), the things will break after some time. Further reading: [Working with Julia projects](https://ufechner7.github.io/2022/08/16/julia-projects.html).

## Which operating systems are supported?
Linux, Windows and Mac are supported. In some of the examples you might have to replace `/` with `\\` on Windows. If you still have a problem, create an issue on [https://github.com/ufechner7/FLORIDyn.jl/issues](https://github.com/ufechner7/FLORIDyn.jl/issues).

## How to create a new simulation?
- in the folder data, copy `2021_9T_Data.yaml` to a new file in the same folder. 
- copy the subfolder `2021_9T_Data` with the `.csv` files to a new folder with the same name as your new `.yaml` file (but without the suffix `.yaml`)
- copy the file `examples/main.jl` to a new file
- edit the line `settings_file = "data/2021_9T_Data.yaml"` to match the name of your new configuration file
- edit the new input files to match your test case
- run your new example using "include("examples/<new_main.jl>")" 

## I have a Problem. Where can I get help?

Ask your question at [Discourse](https://discourse.julialang.org/). Most of the times you will get an answer in 15 min. Half of the people who answer are scientists.