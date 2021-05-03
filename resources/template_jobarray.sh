#!/bin/bash
                                                                                                                                                                                              
#SBATCH -p [SPECIFY PARTITION NAME]                     # queue                                                                                                                          
#SBATCH --time=8:30:00                                  # wall-clock time (mins:secs)                                                                                                           
#SBATCH -c 1                                            # requesting 1 core                                                                                                                    
#SBATCH --mem=4000M                                                                                                                                                                             
#SBATCH -o /dev/null                                    # File to which STDOUT + STDERR will be written 
