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

#include "povwriter.h" 

// globals 

bool config_dom_initialized = false; 
pugi::xml_document config_doc; 	
pugi::xml_node config_root; 

Options options; 
std::vector<Cell_Colors> cell_color_definitions; 

void (*pigment_and_finish_function)(Cell_Colorset&,std::vector<std::vector<double>>&,int); 

std::string VERSION = "1.0.0"; 

void plot_cell( std::ostream& os, std::vector<std::vector<double>>& MAT, int i )
{
	// bookkeeping 
	static Cell_Colorset colors; 
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
	static double nuclear_offset = options.nuclear_offset; 
	
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
	static double bound = options.cell_bound;  
	
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

void cancer_immune_pigment_and_finish_function( Cell_Colorset& colors, std::vector<std::vector<double>>& MAT, int i ) 
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

void standard_pigment_and_finish_function( Cell_Colorset& colors, std::vector<std::vector<double>>& MAT, int i ) 
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

Cell_Colorset::Cell_Colorset( )
{
	cyto_pigment = {1,1,1,0}; 
	nuclear_pigment = {.125,.125,.125,0};
	
	finish = {0.05,1,0.1};
	
	return; 
}

Cell_Colors::Cell_Colors( )
{
	type = 0; 
	return; 
}

void setup_cell_color_definitions( void )
{
	cell_color_definitions.resize( 1 ); 
	
	return; 
}

	
bool load_config_file( std::string filename )
{
	std::cout << "Using config file " << filename << " ... " << std::endl ; 
	pugi::xml_parse_result result = config_doc.load_file( filename.c_str()  );
	
	if( result.status != pugi::xml_parse_status::status_ok )
	{
		std::cout << "Error loading " << filename << "!" << std::endl; 
		return false;
	}
	
	config_root = config_doc.child("povwriter_settings");
	config_dom_initialized = true; 
	
	pugi::xml_node node; 
	
	// figure out filenames, etc. 
	
	node = xml_find_node( config_root , "save" );
	
	options.folder = xml_get_string_value( node, "folder" ) ;
	options.filebase = xml_get_string_value( node, "filebase" ) ;
	options.time_index = xml_get_int_value( node, "time_index" ) ; 
	
	char temp [1024]; 
	sprintf( temp , "./%s/%s%08i_cells_physicell.mat" , options.folder.c_str(), options.filebase.c_str() , options.time_index );
	options.filename = temp; 

	std::cout << "Processing file " << options.filename << "... " << std::endl; 
	
	// decide which function to use 
	
	node = xml_find_node( config_root , "options" ); 
	if( xml_get_bool_value( node , "use_standard_colors" ) == true )
	{
		std::cout << "\tUsing standard coloring function ... "<< std::endl; 
		pigment_and_finish_function = standard_pigment_and_finish_function; 
	}
	else
	{
		std::cout << "\tUsing user-defined coloring in my_pigment_and_finish_function ... " << std::endl; 
		pigment_and_finish_function = my_pigment_and_finish_function; 		
	}
	options.nuclear_offset = xml_get_double_value( node, "nuclear_offset" ); 
	options.cell_bound = xml_get_double_value( node, "cell_bound" ); 
		
	// now, set clipping planes 
	
	node = xml_find_node( config_root , "clipping_planes" ); 
	
	int i = 0; 
	Clipping_Plane cp; 
	
	pugi::xml_node node1 = node.first_child(); 
	while( node1 )
	{
		// get the CSV 
		std::string temp = xml_get_my_string_value( node1 ) ;
		
		// convert to a vector 
		csv_to_vector( temp.c_str() , cp.coefficients ); 
		
		// add the clipping plane 
		default_POV_options.clipping_planes.push_back( cp ); 
		
		// find teh next clipping plane 
		node1 = node1.next_sibling(); 
		i++; 
	}
	std::cout << "Found " << i << " clipping planes ... " << std::endl; 
	
	// read all the cell color definitions 
	
	node = xml_find_node( config_root , "cell_color_definitions" ); 
	
	i = 0; 
	node1 = node.first_child(); 	
	
	while( node1 )
	{
		cell_color_definitions.resize( cell_color_definitions.size()+1 ); 
		
		// set type 
		
		cell_color_definitions[i].type = atoi( node1.attribute( "type" ).value() );
	
		// live 
		node = xml_find_node( node1 , "live" ); 
		std::string temp = xml_get_string_value( node, "cytoplasm" ); 
		csv_to_vector( temp.c_str() , cell_color_definitions[i].live.cyto_pigment ); 
		
		temp = xml_get_string_value( node, "nuclear" ); 
		csv_to_vector( temp.c_str() , cell_color_definitions[i].live.nuclear_pigment ); 

		temp = xml_get_string_value( node, "finish" ); 
		csv_to_vector( temp.c_str() , cell_color_definitions[i].live.finish ); 
		node = node.parent(); 
		
		// apoptotic 
		node = xml_find_node( node1 , "apoptotic" ); 
		temp = xml_get_string_value( node, "cytoplasm" ); 
		csv_to_vector( temp.c_str() , cell_color_definitions[i].apoptotic.cyto_pigment ); 
		
		temp = xml_get_string_value( node, "nuclear" ); 
		csv_to_vector( temp.c_str() , cell_color_definitions[i].apoptotic.nuclear_pigment ); 

		temp = xml_get_string_value( node, "finish" ); 
		csv_to_vector( temp.c_str() , cell_color_definitions[i].apoptotic.finish ); 
		node = node.parent(); 
		
		// necrotic 
		node = xml_find_node( node1 , "necrotic" ); 
		temp = xml_get_string_value( node, "cytoplasm" ); 
		csv_to_vector( temp.c_str() , cell_color_definitions[i].necrotic.cyto_pigment ); 
		
		temp = xml_get_string_value( node, "nuclear" ); 
		csv_to_vector( temp.c_str() , cell_color_definitions[i].necrotic.nuclear_pigment ); 

		temp = xml_get_string_value( node, "finish" ); 
		csv_to_vector( temp.c_str() , cell_color_definitions[i].necrotic.finish ); 
		node = node.parent(); 
		
		node1 = node1.next_sibling(); 
		i++; 
	}
	std::cout << "Found " << cell_color_definitions.size()  << " cell color definitions ... " << std::endl; 
	
	// set up camera 
	
	node = xml_find_node( config_root , "camera" );
	options.camera_distance = xml_get_double_value( node, "distance_from_origin" ); 
	options.camera_theta = xml_get_double_value( node, "xy_angle" ); 
	options.camera_phi = xml_get_double_value( node, "yz_angle" ); 
	
	return true; 	
}

Options::Options()
{
	folder = "sample"; 
	filebase = "output"; 
	time_index = 3696; 
	filename = "./sample/output0003696_physicell_cells.mat" ; 

	double pi = 3.141592653589793;

	camera_distance = 1500; 
	camera_theta = 5*pi/4; 
	camera_phi = pi/3; 
	
	nuclear_offset = 0.1; 
	
	cell_bound = 750; 
	
	return; 
}

void my_pigment_and_finish_function( Cell_Colorset& colors, std::vector<std::vector<double>>& MAT, int i )
{
	cancer_immune_pigment_and_finish_function(colors,MAT,i); 
	return; 
}

void display_splash( std::ostream& os )
{
	os << "povwriter version " << VERSION << std::endl << std::endl; 
	
	os << "Usage: " << std::endl 
	   << "================================================================================" << std::endl 
		<< "povwriter\t\t: " << "run povwriter with config file ./config/settings.xml" << std::endl << std::endl

		<< "povwriter FILENAME.xml\t: " << "run povwriter with config file FILENAME.xml" << std::endl << std::endl

		<< "povwriter x:y:z\t\t: " << "run povwriter on data in FOLDER with indices from x " << std::endl 
		<< "               \t\t  " << "to y in incremenets of z" << std::endl << std::endl 
		<< "               \t\t " << "Example: ./povwriter 0:2:10 processes files: " << std::endl 
		<< "               \t\t " << "         ./FOLDER/FILEBASE00000000_physicell_cells.mat" << std::endl 
		<< "               \t\t " << "         ./FOLDER/FILEBASE00000002_physicell_cells.mat" << std::endl 
		<< "               \t\t " << "         ..." << std::endl 
		<< "               \t\t " << "         ./FOLDER/FILEBASE00000010_physicell_cells.mat" << std::endl 
		<< "               \t\t " << "(See the config file to set FOLDER and FILEBASE)" << std::endl << std::endl 
		
		
		<< std::endl; 

	return; 
}

bool is_xml( std::string filename )
{
	if( strstr( filename.c_str() , ".xml" ) )
	{ return true; }
	return false;
}

bool is_xml( char* filename )
{
	if( strstr( filename, ".xml" ) )
	{ return true; }
	return false;
}

std::vector<int> create_index_list( char* input )
{
	std::vector<int> output; 
	if( strlen( input ) < 1 )
	{ return output; } 
	
	// process list of the form x:y:z 
	char* search = input; 
	char* pch = strstr( search, ":" ); 
	if( pch )
	{
		int start = atoi( search ); 
		int increment = atoi( pch+1 ); 
		pch = strstr( pch+1 , ":" ); 
		int end = atoi( pch+1 ); 
		
		for( int i=start ; i <= end ; i += increment )
		{
			output.push_back( i ); 
		}
		return output; 
	}
	
	// process list in the form of a comma-separated list 
	search = input; 
	pch = strstr( search, "," ); 
	if( pch )
	{
		output.push_back( atoi(input) ); 
		while( pch ) 
		{
			output.push_back( atoi(pch+1) ); 
			pch = strstr( pch+1 , "," ); 
		}
		return output; 
	}
	
	// process list in the form of a space-separated list 
	search = input; 
	pch = strstr( search, " " ); 
	if( pch )
	{
		output.push_back( atoi(input) ); 
		while( pch ) 
		{
			output.push_back( atoi(pch+1) ); 
			pch = strstr( pch+1 , " " ); 
		}
		return output; 
	}	
	
	return output; 
}


