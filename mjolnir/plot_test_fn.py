import subprocess as sub
import mjolnir_plot_helper as mph
import argparse

#script designed to test mjolnir/regrid/hamarr functions

simulations_path = 'simulations/repo_benchmarks'
fig_destination = 'simulations/repo_benchmarks/test_figures'

def moveit(file,prefix,dest):
    #moves plots to new destination
    pname = file.split('/')[-1]
    sub.run(['cp '+file+' '+dest+'/'+prefix+'_'+pname],shell=True)
