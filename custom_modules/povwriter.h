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

#include "../modules/PhysiCell_POV.h"
#include "../modules/PhysiCell_pugixml.h"
#include "../BioFVM/BioFVM_matlab.h" 
#include "../BioFVM/BioFVM_vector.h" 

using namespace BioFVM; 
using namespace PhysiCell; 

extern bool config_dom_initialized; 
extern pugi::xml_document config_doc; 	
extern pugi::xml_node config_root; 

class Options
{
 public:
	std::string folder; 
	std::string filebase; 
	std::string filename; 
	int time_index; 
	
	double camera_distance; 
	double camera_theta;
	double camera_phi; 
	
	double nuclear_offset;
	double cell_bound; 
	
	Options(); 
};

extern Options options; 

class Cell_Colorset
{
 public:
	std::vector<double> cyto_pigment; 
	std::vector<double> nuclear_pigment; 
	std::vector<double> finish; 
	
	Cell_Colorset(); 
}; 

class Cell_Colors
{
 public:
	int type; 
	Cell_Colorset live; 
	Cell_Colorset apoptotic; 
	Cell_Colorset necrotic; 
	
	Cell_Colors(); 
};

extern std::vector<Cell_Colors> cell_color_definitions; 

bool load_config_file( std::string filename ); 
void setup_cell_color_definitions( void ); 

extern void (*pigment_and_finish_function)(Cell_Colorset&,std::vector<std::vector<double>>&,int); 
	
void cancer_immune_pigment_and_finish_function( Cell_Colorset& colors, std::vector<std::vector<double>>& MAT, int i ); 
void standard_pigment_and_finish_function( Cell_Colorset& colors, std::vector<std::vector<double>>& MAT, int i );  
void my_pigment_and_finish_function( Cell_Colorset& colors, std::vector<std::vector<double>>& MAT, int i ); 

void plot_cell( std::ostream& os, std::vector<std::vector<double>>& MAT, int i );

void plot_all_cells( std::ostream& os , std::vector<std::vector<double>>& MAT );

