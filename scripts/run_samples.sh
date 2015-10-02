#!/bin/bash

#simple_run -in zp_atlfast_50ns.txt -out mini_atlfast_50
#simple_run -in zp_fullsim_50ns.txt -out mini_fullsim_50
#simple_run -in zp_atlfast_25ns.txt -out mini_atlfast_25
#simple_run -in zp_fullsim_25ns.txt -out mini_fullsim_25

#simple_run -in zp_atlfast_25ns.txt -out atlfast_25 -atlFast 1
#simple_run -in zp_atlfast_50ns.txt -out atlfast_50 -atlFast 1

simple_run -in zp_fullsim_25ns.txt -out fullsim_25
simple_run -in zp_fullsim_50ns.txt -out fullsim_50
