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
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
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
#include "./BioFVM/BioFVM_matlab.h" 

using namespace BioFVM; 

class cell_colorset
{
 public:
	std::vector<double> cyto_pigment; 
	std::vector<double> nuclear_pigment; 
	std::vector<double> finish; 
	
	cell_colorset(); 
}; 

class cell_colors
{
 public:
	int type; 
	cell_colorset live; 
	cell_colorset apoptotic; 
	cell_colorset necrotic; 
	
	cell_colors(); 
};

std::vector<cell_colors> cell_color_definitions; 

void setup_cell_color_definitions( void ); 

void (*pigment_and_finish_function)(cell_colorset&,std::vector<std::vector<double>>&,int); 
	
void cancer_immune_pigment_and_finish_function( cell_colorset& colors, std::vector<std::vector<double>>& MAT, int i ); 

void standard_pigment_and_finish_function( cell_colorset& colors, std::vector<std::vector<double>>& MAT, int i );  

void plot_cell( std::ostream& os, std::vector<std::vector<double>>& MAT, int i );

void plot_all_cells( std::ostream& os , std::vector<std::vector<double>>& MAT );

int main( int argc, char* argv[] )
{
	char temp [1024]; 
	sprintf( temp , "./cancer_immune_3D/output%08i_cells_physicell.mat" , atoi( argv[1] ) );
	std::cout << "Processing file " << temp << "... " << std::endl; 

	// read the matrix 
	std::string filename = temp; 
	std::vector< std::vector<double> > MAT = read_matlab( filename );
	std::cout << MAT.size() << " x " << MAT[0].size() << std::endl; 
	
	// set options 
	
	pigment_and_finish_function = cancer_immune_pigment_and_finish_function; 
	
	Clipping_Plane cp; 

	cp.coefficients = {0,-1,0,0};
	default_POV_options.clipping_planes.push_back( cp ); 

	cp.coefficients = {-1,0,0,0};
	default_POV_options.clipping_planes.push_back( cp ); 

	cp.coefficients = {0,0,1,0};
	default_POV_options.clipping_planes.push_back( cp ); 
	
	double pi = 3.141592653589793;
	
	default_POV_options.set_camera_from_spherical_location( 1500, 5*pi/4.0 , pi/3.0 ); // do
	default_POV_options.light_position[0] *= 0.5; 

	
	// start output 
	sprintf( temp , "pov%08i.pov" , atoi( argv[1] ) ); 
	std::ofstream os( temp, std::ios::out ); 
	
	
	std::cout << "Creating file " << temp << " for output ... " << std::endl; 
	Write_POV_start( os ); 
	
	std::vector<double> pigment = {1,0,0,0}; 
	std::vector<double> finish = {0.1,0.8,0.5}; 
	
	std::vector<double> center = {0,0,0};
	double radius = 500.0; 
	
	double temp_constant = 0.238732414637843; // 3/(4*pi)
	
	// now, place the cells	
	std::cout << "Writing " << MAT[0].size() << " cells ... " <<std::endl; 
	
	plot_all_cells(os,MAT);
	os.close(); 
	
	std::cout << "done!" << std::endl; 
	
	return 0;
}

void plot_cell( std::ostream& os, std::vector<std::vector<double>>& MAT, int i )
{
	// bookkeeping 
	static cell_colorset colors; 
	static bool setup_done = false; 
	if( setup_done == false )
	{
		colors.cyto_pigment = {1,1,1,0}; 
		colors.nuclear_pigment = {.5,.5,.5,0};
		colors.finish = {0.05,1,0.1};
		setup_done = true; 
	}
	
	static std::vector<double> center = {0,0,0};
	static double radius; 
	static double temp_constant = 0.238732414637843; // 3/(4*pi)
	
	// get position 
	
	center[0] = MAT[1][i]; 
	center[1] = MAT[2][i]; 
	center[2] = MAT[3][i]; 
	
	// first, plot the cytoplasm 
		
	radius = pow( temp_constant * MAT[4][i] , 0.33333333333333333333333333333 ); 
		
	// test against clipping planes 
	
	bool render = true; 
	if( default_POV_options.clipping_planes.size() > 0 )
	{ render = false; }

	std::vector<int> intersection_indices = {}; 

	double testval; 
	bool intersect = false; 
	
	for( int i=0; i < default_POV_options.clipping_planes.size() ; i++ )
	{
		testval = default_POV_options.clipping_planes[i].signed_distance_to_plane( center ); 
		
		if( testval <= -radius )
		{ render = true; }
		if( testval > -radius && testval <= radius )
		{
			render = true; 
			intersect = true; 
			intersection_indices.push_back( i ); 
		}
	}
	
	if( intersect )
	{ os << "intersection{ " << std::endl ; }
	
	if( render )
	{
		// alt_pigment_and_finish_function( cyto_pigment, nuclear_pigment, finish, MAT, i ); 
		pigment_and_finish_function( colors, MAT, i ); 

		if( intersect )
		{
			// if( intersection_indices.size() > 1 )
			{ os << "union{ " << std::endl ; }
			
			int i; 
			for( int i=0; i < default_POV_options.clipping_planes.size() ; i++ )
			{
				os	<< "plane{<" << default_POV_options.clipping_planes[i].coefficients[0] << "," 
					<< default_POV_options.clipping_planes[i].coefficients[1] << "," 
					<< default_POV_options.clipping_planes[i].coefficients[2] << ">, " 
					<< default_POV_options.clipping_planes[i].coefficients[3] << std::endl 
					<< " pigment {color rgbf<" 
						<< colors.cyto_pigment[0] << "," 
						<< colors.cyto_pigment[1] << "," 
						<< colors.cyto_pigment[2] << "," 
						<< colors.cyto_pigment[3] << ">}" << std::endl
					<< " finish {ambient " << colors.finish[0] 
					<< " diffuse " << colors.finish[1] 
					<< " specular " << colors.finish[2] << "} }" << std::endl;
			}
			
			// if( intersection_indices.size() > 1 )
			{ os << "}" << std::endl; }
		}
		
		Write_POV_sphere( os, center, radius, colors.cyto_pigment, colors.finish ); 
	}

	if( intersect )
	{ os << "}" << std::endl; }

	// now, plot the nucleus 
	
	radius = pow( temp_constant * MAT[9][i] , 0.33333333333333333333333333333 ); 
		
	// test against clipping planes 
	
	render = true; 
	if( default_POV_options.clipping_planes.size() > 0 )
	{ render = false; }

	intersection_indices.resize(0); 

	intersect = false; 
	
	// offset the nuclear clipping just tiny bit, to avoid 
	// graphical artifacts where the cytoplasm and nucleus 
	// blend into each other 
	double nuclear_offset = 0.1; 
	
	for( int i=0; i < default_POV_options.clipping_planes.size() ; i++ )
	{
		testval = default_POV_options.clipping_planes[i].signed_distance_to_plane( center ); 
		
		if( testval <= -(radius+nuclear_offset) )
		{ render = true; }
		if( testval > -(radius+nuclear_offset) && testval <= (radius+nuclear_offset) )
		{
			render = true; 
			intersect = true; 
			intersection_indices.push_back( i ); 
		}
	}
	
	if( intersect )
	{ os << "intersection{ " << std::endl ; }
	
	if( render )
	{

		if( intersect )
		{
			// if( intersection_indices.size() > 1 )
			{ os << "union{ " << std::endl ; }
			
			int i; 
			for( int i=0; i < default_POV_options.clipping_planes.size() ; i++ )
			{
				os	<< "plane{<" << default_POV_options.clipping_planes[i].coefficients[0] << "," 
					<< default_POV_options.clipping_planes[i].coefficients[1] << "," 
					<< default_POV_options.clipping_planes[i].coefficients[2] << ">, " 
					<< default_POV_options.clipping_planes[i].coefficients[3]+nuclear_offset << std::endl 
					<< " pigment {color rgbf<" 
						<< colors.nuclear_pigment[0] << "," 
						<< colors.nuclear_pigment[1] << "," 
						<< colors.nuclear_pigment[2] << "," 
						<< colors.nuclear_pigment[3] << ">}" << std::endl
					<< " finish {ambient " << colors.finish[0] 
					<< " diffuse " << colors.finish[1] 
					<< " specular " << colors.finish[2] << "} }" << std::endl;
			}
			
			// if( intersection_indices.size() > 1 )
			{ os << "}" << std::endl; }
		}
		default_POV_options.no_shadow = true; 
		Write_POV_sphere( os, center, radius, colors.nuclear_pigment, colors.finish ); 
		default_POV_options.no_shadow = false; 
	}

	if( intersect )
	{ os << "}" << std::endl; }

	return; 
}

void plot_all_cells( std::ostream& os , std::vector<std::vector<double>>& MAT )
{
	double bound = 750.0; 
	
	for( int i = 0 ; i < MAT[0].size() ; i++ )
	{
		if( MAT[1][i] > -bound && MAT[1][i] < bound &&
		MAT[2][i] > -bound && MAT[2][i] < bound &&
		MAT[3][i] > -bound && MAT[3][i] < bound )
		{		
			plot_cell( os, MAT, i ); 
		}
	}	

	return; 
}

void cancer_immune_pigment_and_finish_function( cell_colorset& colors, std::vector<std::vector<double>>& MAT, int i ) 
{
	// first, some housekeeping
	static int type_index = 5; //column that stores cell type (integer)
	static int cycle_model_index = 6; // column that stores which cycle (or death) model (integer) 
	static int data_index = 27; // 
	
	colors.finish = { 0.025 , 1 , 0.1 }; 
	
	// if this is an immune cell, make it red
	if( (int) MAT[type_index][i] == 1 )
	{
		colors.cyto_pigment = {1.0, 0.0, 0.0 , 0.0};  
		colors.nuclear_pigment = {0,0.125,0.0 , 0.0};  
		
		return; 
	}
	
	bool necrotic = false; 
	bool apoptotic = false; 
	bool live = true; 
	int cycle_model = (int) round( MAT[cycle_model_index][i] ); 
	if( cycle_model == 100 )
	{
		apoptotic = true;
		live = false;
	}		
	if( cycle_model == 101 )
	{
		necrotic = true;
		live = false;
	}		

	// live cells are green, but shaded by oncoprotein value 
	if( live == true )
	{
		double oncoprotein = MAT[data_index][i]; // 0.5 * MAT[27][i];  
		
		// map [0.5 1.5] to [0 1]
		if( oncoprotein > 1.5 )
		{ oncoprotein = 1.5; }
		if( oncoprotein < 0.5 )
		{ oncoprotein = 0.5; } 
		oncoprotein -= 0.5; 
		
		colors.cyto_pigment[0] = oncoprotein;
		colors.cyto_pigment[1] = oncoprotein; 
		colors.cyto_pigment[2] = 1.0 - oncoprotein;
		
		colors.nuclear_pigment[0] = colors.cyto_pigment[0] * 0.125; 
		colors.nuclear_pigment[1] = colors.cyto_pigment[1] * 0.125; 
		colors.nuclear_pigment[2] = colors.cyto_pigment[2] * 0.125; 
		
		//	sprintf( szTempString , "rgb(%u,%u,%u)", (int)round(output[0][0]/2.0) , (int)round(output[0][1]/2.0) , (int)round(output[0][2]/2.0) );
		return; 
	}	
	
	// apoptotic cells
	
	if( apoptotic )
	{
		colors.cyto_pigment = {0,1,1 , 0.0};  
		colors.nuclear_pigment = {0,0.125,0.125 , 0.0};  
	
		return; 
	}
	
	// necrotic cells 
	
	if( necrotic )
	{
		colors.finish = {0.01,0.5,0.1}; 
		colors.cyto_pigment = {1,0.5412,0.1490, 0.0}; 
		colors.nuclear_pigment = {0.125,0.06765,0.018625, 0.0};  
		
		return; 
	}
	
	return; 
}

void standard_pigment_and_finish_function( cell_colorset& colors, std::vector<std::vector<double>>& MAT, int i ) 
{
	// first, some housekeeping
	static int type_index = 5; //column that stores cell type (integer)
	static int cycle_model_index = 6; // column that stores which cycle (or death) model (integer) 
	
	// Search the (global) array of colors 
	// for the cell's type. If it's not found,
	// default to 0. 
	
	int color_index = 0; 
	int cell_type = (int) MAT[type_index][i]; 
	for( int j=0; j < cell_color_definitions.size(); j++ )
	{
		if( cell_type == cell_color_definitions[j].type )
		{ color_index = j; }
	}
	
	// next, see if it's live, apoptotic, or necrotic 
	
	bool necrotic = false; 
	bool apoptotic = false; 
	bool live = true; 
	int cycle_model = (int) round( MAT[cycle_model_index][i] ); 
	
	if( cycle_model == 100 )
	{
		apoptotic = true;
		live = false;
	}		
	if( cycle_model == 101 )
	{
		necrotic = true;
		live = false;
	}		

	// cell is live. use the appropriate colors 
	if( live == true )
	{
		colors = cell_color_definitions[color_index].live; 
		return; 
	}
	
	if( apoptotic == true )
	{
		colors= cell_color_definitions[color_index].apoptotic; 
		return; 
	}
		
	if( necrotic == true )
	{
		colors= cell_color_definitions[color_index].necrotic; 
		return; 
	}
		
	return; 
}

cell_colorset::cell_colorset( )
{
	cyto_pigment = {1,1,1,0}; 
	nuclear_pigment = {.125,.125,.125,0};
	
	finish = {0.05,1,0.1};
	
	return; 
}

cell_colors::cell_colors( )
{
	type = 0; 
	return; 
}

void setup_cell_color_definitions( void )
{
	cell_color_definitions.resize( 1 ); 
	
	return; 
}


