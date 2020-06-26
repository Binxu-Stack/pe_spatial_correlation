# variables to build initial simulation box
# info from website http://www.periodictable.com/Elements/040/data.html

variable Nevery  equal 100
variable Nrepeat  equal 200
variable Nfreq  equal 20000

variable max_r    equal 10.0  # maximum radius for sc


variable T     equal 1400

#variable seed                   equal ${my_seed}
variable seed                   equal 1
boundary         p p p
units            metal
dimension        3
atom_style       atomic

read_restart  ./T_${T}.restart


pair_style  eam/alloy
pair_coeff  * * ZrCu.lammps.eam Zr Cu


# ensamble (NVT)
variable tstep equal "2e-3"
variable Pdamp equal "v_tstep * 1000"  # suggusted value for Nose-Hoover thermostate
variable Tdamp equal "v_tstep * 100"   # suggusted value for Nose-Hoover thermostate

reset_timestep   0
timestep         ${tstep}

#                     1   2       3     4      5   6   7   8   9   10 11  12  13
thermo_style custom step temp ke etotal vol  pe pxx pyy pzz pxy pyz pxz lx
thermo 100

neighbor 2.0 bin
neigh_modify every 4 delay 0 check yes
# give new velocity
#velocity all create ${temp0} ${seed} mom yes rot yes dist gaussian units box

variable comm_cutoff equal  ${max_r}+2  # suppose neighbor skin <= 2.0
comm_modify cutoff ${comm_cutoff}   # ghost atom commucation radius  should be larger than (max_r + neighbour skin)

#compute PE all pe/atom
compute SC all sc 100 1 1 1 2 2 2 cutoff ${max_r}
fix AVE all ave/time ${Nevery} ${Nrepeat} ${Nfreq} c_SC[*] file sc_all.dat mode vector


fix NVT all nvt temp ${T} ${T} ${Tdamp}


run ${Nfreq}
