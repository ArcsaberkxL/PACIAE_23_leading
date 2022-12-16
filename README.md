<!-- This is a README file for usage of PACIAE.
     Written by Markdown language.
                 By Anke at CCNU on 10/16/2022 
 -->


# The parton and hadron cascade model PACIAE

The PACIAE model is a multipurpose Monte Carlo event generator developed to describe a wide range of collisions, including hadron-hadron interactions, hadron interactions off nuclei, and nucleus-nucleus collisions. It is built based on PYTHIA-6.4 and incorporates parton and hadron rescattering stages to take care of the nuclear medium effects.

## Install
Just install the file and decompress it. Then you can get the PACIAE source code directly.


## Usage

We encourage users to run the program on LINUX.
There are two ways to run the program.
 1. Use the normal compiler to complie the soruce code and run the PACIAE program. Take gfortran on LINUX as an example:
    - Compile and link the programs by the command:
        ```
            gfortran -O -C *.f -o paciae.x
        ```
    - Modify the input file of usu.dat according to your wish.
    - Run the program by the command:
        ```
            time ./paciae.x
        ```
      One could use the following command to run the program to in the background and record the log information.
        ```
            nohup time ./paciae.x > paciae.log &
        ```
 2. Use the PACIAE.sh shell-script to compile the code, generate usu.dat file and run program automatically.
    - Modify the PACIAE.sh file as needed.
    - Give executable permission to the PACIAE.sh file by command:
        ```
            chmod +x PACIAE.sh
        ```
    - Run the PACIAE.sh script by the command:
        ```
            time ./PACIAE.sh
        ```
         The more recommended command is:
         ```
            time ./PACIAE.sh | tee $(date "+%Y%m%d%H%M%S").log
         ``` 
         which records the screen information to a log file with date and time.
         Of course, one could run it in the background:
         ```
            nohup time ./PACIAE.sh | tee $(date "+%Y%m%d%H%M%S").log &
         ```
    It is also worth mentioning that tasks can be submitted to computer clusters and supercomputers using the PACIAE.sh script. (It is currently only available in SLRUM and LSF scheduling systems) More detailed information and usage please read the PACIAE.sh file.


## Maintainers

[@ArcsaberHep](https://github.com/ArcsaberHep).


## Contributing

Feel free to dive in! Any bug reports, comments and suggestions are welcome. Please do not hesitate to contact us.


## Contributors

 - [An-Ke Lei](https://inspirehep.net/authors/1965068), ankeleihep@gmail.com or ankelei@mains.ccnu.edu.cn <!-- Key Laboratory of Quark and Lepton Physics (MOE) and Institute of Particle Physics, Central China Normal University, Wuhan 430079, China. --> 
 - [Dai-Mei Zhou](https://inspirehep.net/authors/1030208), zhoudm@mail.ccnu.edu.cn <!-- Key Laboratory of Quark and Lepton Physics (MOE) and Institute of Particle Physics, Central China Normal University, Wuhan 430079, China. --> 
 - [Yu-Liang Yan](https://inspirehep.net/authors/1051028), yanyl@ciae.ac.cn <!-- China Institute of Atomic Energy, P.O. Box 275 (10), Beijing, 102413,China. -->
 - [Ben-Hao Sa](https://inspirehep.net/authors/990834), sabh@ciae.ac.cn (model founder) <!-- China Institute of Atomic Energy, P.O. Box 275 (10), Beijing, 102413,China. -->


## License

[GPL v2.0](LICENSE)


## Released papers

 - Revisiting the centrality definition and observable centrality dependence of relativistic heavy-ion collisions in PACIAE model, [Comput. Phys. Commun. 284 (2023) 108615](https://doi.org/10.1016/j.cpc.2022.108615) or [arXiv:2212.04087 [nucl-th]](https://doi.org/10.48550/arXiv.2212.04087)

 - PACIAE 2.2.1: An updated issue of the parton and hadron cascade model PACIAE 2.2, [Comput. Phys. Commun. 274 (2022) 108289](https://doi.org/10.1016/j.cpc.2022.108289) <!-- or [arXiv: [nucl-th]](https://doi.org/) -->

 - Announcement for the replacement of the PACIAE 2.1 and PACIAE 2.2 series, [Comput. Phys. Commun. 224 (2018) 417-418](https://doi.org/10.1016/j.cpc.2017.10.006) <!-- or [arXiv: [nucl-th]](https://doi.org/) -->

 - An upgraded issue of the parton and hadron cascade model, PACIAE 2.2, [Comput. Phys. Commun. 193 (2015) 89-94](https://doi.org/10.1016/j.cpc.2015.01.022) or [arXiv:1412.7579 [nucl-th]](https://doi.org/10.48550/arXiv.1412.7579)

 - PACIAE 2.1: An updated issue of the parton and hadron cascade model PACIAE 2.0, [Comput. Phys. Commun. 184 (2013) 1476-1479](https://doi.org/10.1016/j.cpc.2012.12.026) or [arXiv:1206.4795 [nucl-th]](https://doi.org/10.48550/arXiv.1206.4795)

 - PACIAE 2.0: An updated parton and hadron cascade model (program) for the relativistic nuclear collisions, [Comput. Phys. Commun. 183 (2012), 333-346](https://doi.org/10.1016/j.cpc.2011.08.021) or [arXiv:1104.1238 [nucl-th]](https://doi.org/10.48550/arXiv.1104.1238)
<br/>

 - Special version: HYDRO-PACIAE, a hydrodynamic and transport hybrid model for ultra-relativistic heavy ion collisions, [J.Phys.G 40 (2013) 025102](https://doi.org/10.1088/0954-3899/40/2/025102) or [arXiv:1110.6704 [nucl-th]](https://doi.org/10.48550/arXiv.1110.6704)


## Update Notes:

<!----------------------------------------------------------------------------->
...(waiting for updating)

12/2022 In version PACIAE2.3

- In subroutine "main" in "main_23.f", the real-time clock randonm seed is introduced.
- In subroutine "scat" in the "parini_23.f", the long-written statement about executing the binary collision by calling PYEVNW / PYEVNT is replaced by a new subroutine "xevent".
- In main_23.f<!-- , eps09.f -->, the statements with old sytanx are re-wrote by more modern-style one.
- The CME is introduced in PACIAE2.3 based on Zhi-Lei's improvement in PACIAE2.2.1b and PACIAE2.2.1c.

<!----------------------------------------------------------------------------->
04/2021 In version PACIAE2.2.1b and PACIAE2.2.1c

- The CME is introduced by Zhi-Lei She et al.

<!-----------------------------------------------------------------------------
01/2020 In version PACIAE2.3

- The eps09 nuclear shadowing is introduced by Liang, whose sebroutine is called "shanul_eps09". An exrta file named "eps09.f" is added.

----------------------------------------------------------------------------->