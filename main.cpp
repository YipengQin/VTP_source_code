#include "stdafx.h"
#include "geodesic_mesh.h"
#include "geodesic_algorithm_exact.h" 

using namespace std;

int main(int argc, char **argv)
{
	if(argc == 1)
	{
		cout <<endl << "OPTIONS: " << endl;
		cout << endl;

		cout << "-m [meshFile]: input model file." << endl << endl;
		cout << "-s [src]: index of source." << endl << endl;
		cout << "-o [outputFile]: output model file." << endl << endl;
		return -1;
	}

	char file_name[255] = {'\0'};
	char outputMesh[255] = {'\0'};
	unsigned source_vertex_index = 0; // Source Vertex

	for (int i = 1; i < argc;)
	{
		if (strcmp(argv[i], "-m") == 0)
		{
			strcpy_s(file_name, argv[i+1]); i+=2;
		}
		else if (strcmp(argv[i], "-s") == 0)
		{
			while (++i < argc && argv[i][0] != '-')
			{
				source_vertex_index = atoi(argv[i]);
				break;
			}
		}
		else if (strcmp(argv[i], "-o") == 0)
		{
			strcpy_s(outputMesh, argv[i+1]); i+=2;
		}
		else ++i;
	}

	std::vector<double> points;	
	std::vector<unsigned> faces;
	std::vector<int> realIndex;
	int originalVertNum = 0;
	
	std::cout << file_name << std::endl;

	// Load Mesh
	bool success = geodesic::read_mesh_from_file(file_name, points, faces, realIndex, originalVertNum);
	if(!success)
	{
		std::cout << "something is wrong with the input file" << std::endl;
		return -2;
	}

	std::cout << "Load Mesh Success..." << std::endl;
	
	// Build Mesh
	geodesic::Mesh mesh;
	mesh.initialize_mesh_data(points, faces);		//create internal mesh data structure including edges

	std::cout << "Build Mesh Success..." << std::endl;
	
	geodesic::GeodesicAlgorithmExact algorithm(&mesh);

	// Propagation
	algorithm.propagate(source_vertex_index);	//cover the whole mesh

	// Print Statistics
	std::cout << endl;
	algorithm.print_statistics();

	// Output Geodesic Distances
	if (outputMesh[0] != '\0')
	{
		ofstream output(outputMesh);

		for(unsigned i=0; i<mesh.vertices().size(); ++i)
		{
			double distance = mesh.vertices()[i].geodesic_distance();

			output << setprecision(20)<< distance <<endl;		//print geodesic distance for every vertex
		}

		output << endl;
		
		output.close();
	}

	std::cout << std::endl;

	return 0;
}