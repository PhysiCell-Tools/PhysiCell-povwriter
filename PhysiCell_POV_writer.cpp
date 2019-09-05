/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2019, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>
#include <string>

#include "./modules/PhysiCell_POV.h"
#include "./modules/PhysiCell_pugixml.h"
#include "./BioFVM/BioFVM_matlab.h" 
#include "./BioFVM/BioFVM_vector.h" 

#include "./custom_modules/povwriter.h" 

int main( int argc, char* argv[] )
{
	display_splash( std::cout ); 
	
	char* my_text = (char*) "0:1:10"; 
	
	std::vector<int> indices = create_index_list( my_text ); 
	std::cout << "["; 
	for( int i=0; i < indices.size() ; i++ )
	{
		std::cout << indices[i] << " " ; 
	}
	
	std::cout << "]" <<std::endl; 
	
	exit(0); 
	
	std::string config_file = "./config/settings.xml"; 
	
	// load and parse settings file(s)
	
	bool XML_status = false; 
	if( argc > 1 )
	{
		if( is_xml(argv[1]) )
		{
			config_file = argv[1]; 
		}
		else
		{
			std::cout << "not an xml file extension" << std::endl; 			
		}
	}
	
	XML_status = load_config_file( config_file ); 
	if( !XML_status )
	{ exit(-1); }
	
	// read the matrix 
	std::vector< std::vector<double> > MAT = read_matlab( options.filename.c_str() );
	std::cout << "Matrix size: " << MAT.size() << " x " << MAT[0].size() << std::endl; 
	
	// set options 
	
	default_POV_options.set_camera_from_spherical_location( options.camera_distance , options.camera_theta, options.camera_phi ); //  1500, 5*pi/4.0 , pi/3.0 ); // do
	default_POV_options.light_position[0] *= 0.5; 
	
	// start output 
	char temp [1024]; 
	sprintf( temp , "pov%08i.pov" , options.time_index ); 
	options.filename = temp ; 
	std::ofstream os( options.filename.c_str() , std::ios::out ); 
	
	std::cout << "Creating file " << options.filename << " for output ... " << std::endl; 
	Write_POV_start( os ); 
	
	// now, place the cells	
	std::cout << "Writing " << MAT[0].size() << " cells ... " <<std::endl; 
	
	plot_all_cells(os,MAT);
	os.close(); 
	
	std::cout << "done!" << std::endl; 
	
	return 0;
}

