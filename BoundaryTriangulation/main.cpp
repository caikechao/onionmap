#include <stdio.h>
#include <string.h>
#include "sensor.h"
#include "common.h"
#include "random.h"
#include "geometry3D.h"


int statistics(double *array, double x);
int QUASI_UBG = 0;
int Number_of_Landmark = 0;
double QUASI_UBG_ALPHA = 0.0;
int LENGTH = 15;

Field *globalField = NULL;

//Ellipsoid
double Low_Edge_Threshold = 0.0;
double T_NS = 54;
double High_Edge_Threshold = 0.0;
//Convex judgment threshold
double Convex_threshold = -0.1;
double Concave_threshold = 0.4;//important
double Concave_threshold1 = 0.5;
double Concave_threshold2 = 0.7;
double Saddle_threshold = 0.15;
int Contain_convex_threshold = 10;
double RATIO = 0;
double COMM_RANGE = 8;
double Gridinterval = 5;
//int KHOP = 4;			//KHOP must be greater than 1
int Total_Number_of_Triangles = 0;
int RemainingEdges = 0;
int UsedLandmarks = 0;
FILE *g_fp;

//Minimum landmark distance 
const int kHOP = 4-1 + 1;

int main(int argc, char** argv)
{
	FILE *fp;
	errno_t err;
	err = fopen_s(&fp, "triangle_vicility", "w");
	fclose(fp);

	err = fopen_s(&fp, "Line", "w");
	fclose(fp);

	err = fopen_s(&fp, "treeRoutRealandShortestPath", "w");
	fclose(fp);

	err = fopen_s(&fp, "SphereRoutRealandShortestPath", "w");
	fclose(fp);

	err = fopen_s(&fp, "triangle", "w");
	fclose(fp);
	globalField = new Field();
	globalField->initDeployment("Initsensors.txt", "TopologyVicility.txt");

	// Calculate sensor node importance.
	//cerr<<"Calculating sensor node importance..."<<endl;
	//globalField->calSensorImportance();
	//cerr<<"Boundary sensor nodes' importances are defined."<<endl;

	// Landmark selection in the boundary nodes.
	cerr<<"Selecting landmark..."<<endl;
	globalField->STARTLEVEL = 1;

	globalField->banLevel = 16;
	globalField->MAX_BOUNDARY_TYPE = 17;
	globalField->PickingLandmarkWithSignificance(kHOP);
	/*cerr<<globalField->LandmarkSet.size()<<" landmarks are selected."<<endl;*/	

	// Generate voronoi cell
	cerr<<"Building voronoi cell..."<<endl;
	globalField->GenerateVoronoiCell();
	cerr<<"Voronoi cell constructed."<<endl;

	// Dump landmarks.
	// 	cerr<<"Dumping landmarks..."<<endl;
	// 	globalField->dumpTopology("topology.landmarks");
	// 	cerr<<"Dumping over."<<endl;

	// Triangulation of the landmarks.
	cerr<<"Finding edges of landmarks..."<<endl;
	globalField->egdesOfLandmarks();
	cerr<<"All legal edges are found."<<endl;

	cerr<< "Start triangulation..."<<endl;
	globalField->landmarkTriangulation();
	cerr<< "Topology is triangulated."<<endl;

	cerr<<"Out put the layers of the triangle..."<<endl;
	globalField->dumpLayersOfTriangle("Layers_Triangles");
	cerr<<"The layers are dumped."<<endl;
	// 	Testing triangle...
	// 	cerr<<"testing output triangle\n";
	// 	globalField->testingTriangleOutPut();
	// 	cerr<<"Triangle file generated.\n";

	// Sphere ricci flow..
	cerr<<"Starting ricci flow..."<<endl;
	globalField->iTypeRicciFlow();
	cerr<<"Flow over."<<endl;

	// 	globalField->dumpTopology("topology");
	cerr<<"Dump ricci location of the landmarks..."<<endl;
	globalField->dumpItypeRicciLocation("topology.ricci");
	cerr<<"Ricci location dumped."<<endl;


	cerr<<"Start spherical mapping..."<<endl;
	globalField->ItypeStereographicProjection();
	cerr<<"Spherical mapping complete."<<endl;


	//the sphere alignment.
// 	cerr<<"Sphere alignment..."<<endl;
// 	globalField->sphereAlignment();
// 	cerr<<"alignment complete.."<<endl;

	// Give a Landmark uplevel and lowlevel neighbors.


	cerr<<"Dump sphere location of the landmarks..."<<endl;
	globalField->dumpItypeSphereLocation("topology.sphere");
	cerr<<"sphere location dumped."<<endl;


	// Start routing...
	cerr<<"Start routing... "<<endl;
	//building the spanning tree.
 	Sensor* root = globalField->sensorPool + globalField->ROOT;
 	cerr<<"The root sensor generate spanning tree..."<<endl;
	root->generateSpanningtree();
	cerr<<"Naming the sensor..."<<endl;
	root->namingSensors();
	globalField->findAChildForSensors();
	cerr<<"Routing testing..."<<endl;
	globalField->CORE_RADIUS = 3;
	int i;
	int nRounds = 5000;
	int actualRounds = 0;

	bool treeRouting = true;
	bool sphereRouting = true;

	double totalTree = 0;
	double totalSphere = 0;

	double maxTree = 0;
	double maxSphere = 0;

	vector<double> stretchTree;
	stretchTree.reserve(nRounds);
	vector<double> stretchSphere;
	stretchSphere.reserve(nRounds);

	double percentTree[6] = {0};
	double percentSphere[6] = {0};
	UniformRandVar randv(0, globalField->nSensors - 1);

	for(int loop = 0; loop < nRounds; loop++)
	{
		int src = (int)randv.value();
		int dest = (int)randv.value();
		Sensor* srcSensor  = globalField->sensorPool + src;
		Sensor* destSensor = globalField->sensorPool + dest;
// 		auto it_src_maxLevel  = max_element(srcSensor->levelsBelongTo.begin(),srcSensor->levelsBelongTo.end());
// 		auto it_dest_maxLevel = max_element(destSensor->levelsBelongTo.begin(),destSensor->levelsBelongTo.end());
// 		int src_maxLevel  = *it_src_maxLevel;
// 		int dest_maxLevel = *it_dest_maxLevel;
		if(src == dest || srcSensor->neighbors.find(destSensor) != srcSensor->neighbors.end()
			/*|| srcSensor->inCommonType(destSensor) > 2 || src_maxLevel >2 || dest_maxLevel >2*/ )
			continue;

		actualRounds++;

		if (loop % 100 == 0)
		{
			cout<<"*";
		}
		double ratioTree = 0;
		double ratioSphere = 0;

		if(treeRouting)
		{

			ratioTree = globalField->tryTreeRoute(src,dest);
			totalTree += ratioTree;

			if (ratioTree > maxTree)
			{
				maxTree = ratioTree;
			}

			stretchTree.push_back(ratioTree);
		}

		if (sphereRouting)
		{

			ratioSphere = globalField->trySphereRoute(src,dest,true);
			totalSphere += ratioSphere;

			if (ratioSphere > maxSphere)
			{
				maxSphere = ratioSphere;
			}

			stretchSphere.push_back(ratioSphere);

		}
		

		statistics(percentTree, ratioTree);
		statistics(percentSphere, ratioSphere);

	}
	cerr<<endl;
	// Print the information.
	for(i = 0; i < 6; i++)
	{
		fprintf(stderr, "\t%5.2f\t%5.2f\t\n", percentTree[i] / actualRounds, percentSphere[i] / actualRounds);
	}
	fprintf(stderr, "avg:\t%5.2f\t%5.2f\t\n", totalTree / actualRounds, totalSphere / actualRounds);
	fprintf(stderr, "max:\t%5.2f\t%5.2f\t\n", maxTree,maxSphere);
	if (treeRouting)
	{
		printf("Tree routing...\n");

		sort(stretchTree.begin(), stretchTree.end());
		int lowidx = int(0.05 * stretchTree.size());
		int upidx = int(0.95 * stretchTree.size());

		printf("Tree routing: 5-th percentile: %lf, 95-th percentile: %lf\n", stretchTree[lowidx], stretchTree[upidx]);
	}
	if (sphereRouting)
	{
		printf("Sphere routing...\n");

		sort(stretchSphere.begin(), stretchSphere.end());
		int lowidx = int(0.05 * stretchSphere.size());
		int upidx = int(0.95 * stretchSphere.size());

		printf("Sphere routing: 5-th percentile: %lf, 95-th percentile: %lf\n", stretchSphere[lowidx], stretchSphere[upidx]);
	}

	if (globalField)
	{
		delete globalField;
	}
	
	char buf[128];
	gets_s(buf);
	return 0;

}

int statistics(double *array, double x)
{

	if (x < 1.1)
	{
		array[0]++;
	}
	else if (x < 1.5)
	{
		array[1]++;
	}
	else if (x < 2.0)
	{
		array[2]++;
	}
	else if (x < 3.0)
	{
		array[3]++;
	}
	else if (x < 5.0)
	{
		array[4]++;
	}
	else
	{
		array[5]++;
	}

	return 0;
}