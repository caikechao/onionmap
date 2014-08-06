#include "sensor.h"
#include "random.h"
#include "assert.h"
#include "common.h"
#include "geometry3D.h"
#include <math.h>
#include <string>
#include <complex>
Field::Field()
{
	nSensors = 0;
	comRadius = 0;
	avg1HopNeighborhoodSize = 0;
	avgNHopNeighborhoodSize = 0;
	avgOverallNodeDegree = 0;
	avgOverallHopDistWithinCircle = 0;
	Avg_distance_err = 0;
	num_of_outer_boundary_surface = 0;
	NumOfBoundaryPart = 0;
	EdgeSet.clear();
	RefineSensorSet.clear();
	VoronoiCellSensorSet.clear();

	removedTris.clear();
	noCrossEdgePool.clear();
	removedSurfaceNodes.clear();
	ROOT = 0;

}

Field::Field(Field &field)
{
	nSensors = field.nSensors;
	comRadius = field.comRadius;
	avg1HopNeighborhoodSize = field.avg1HopNeighborhoodSize;
	avgNHopNeighborhoodSize = field.avgNHopNeighborhoodSize;
	avgOverallNodeDegree = field.avgOverallNodeDegree;
	sensorPool = new Sensor[nSensors];
	for(int i = 0; i < nSensors; ++i)
	{
		sensorPool[i] = field.sensorPool[i];
	}
	EdgeSet.assign(field.EdgeSet.begin(), field.EdgeSet.end());
}

Field::~Field()
{
	delete [] sensorPool;

}

//read file
void Field::initDeployment(const char* topoFile1, const char* topoFile2)
{
	int i = 0;
	vector<Point> validPoints;
	
	if(topoFile1) 
	{
		FILE *fp;
		errno_t err;
		err = fopen_s(&fp, topoFile1, "r");
		//int debugN = 0;

		char line[128];
		if (err)  
		{
			fprintf_s(stderr, "no topo file %s found\n", topoFile1);
			exit(0);
		}

		Point s;
		int countIndex = 0;
		do 
		{
			
			memset(line, 0, sizeof(line));
			fgets(line, sizeof(line), fp);

			if((*line) == '#')
			{
				sscanf_s(line, "#%*d %lf %d %d", &COMM_RANGE, &MAX_BOUNDARY_TYPE,&ROOT);
				cout<<"Root is "<<ROOT<<endl;
				printf("Communication range is %lf\nMax level is %d\n", COMM_RANGE, MAX_BOUNDARY_TYPE);
			}
			else if ((*line) == '~')// Read the levels 
			{
 			//	cout<<debugN++<<endl;
				s.levelTypes.clear();
				int i = 1;// From the second character...not '~'
				int ilevel;
				
				char val[1024];
			    memset(val, 0, sizeof(val));
				while(line[i] != '\n')
				{
					int j = 0;
					while(line[i] != ' ' && line[i] != '\n')
					{
						val[j] = line[i];
						++i;
						++j;
					}
					ilevel = atoi(val);
					memset(val, 0, sizeof(val));
					
					s.levelTypes.insert(ilevel);
					++i;
				}

				validPoints.push_back(s);

			}
			else if(!isspace(*line) && (*line)) 
			{
				sscanf_s(line, "%*d %lf %lf %lf", &s.x, &s.y, &s.z);
			}
		}while(!feof(fp));
		fclose(fp);
		nSensors = validPoints.size();
		comRadius = COMM_RANGE;
		printf("nSensor = %d\n", nSensors);
	}

	sensorPool = new Sensor[nSensors];
	cerr<<"Reading the sensors..."<<endl;
	for(i = 0; i < (signed)validPoints.size(); ++i)
	{
		//initialize index, location, boundaryType
		sensorPool[i].init(i, validPoints[i]);
		sensorPool[i].field = this;

	}
	cerr<<"Reading the neighbors..."<<endl;
	if(topoFile2) 
	{
		FILE *fp;
		errno_t err;
		err = fopen_s(&fp, topoFile2, "r");
		if (err)  
		{
			fprintf_s(stderr, "no topo file %s found\n", topoFile2);
			exit(0);
		}
		char line[1024];
		int count = 0;
		do 
		{
			memset(line, 0, sizeof(line));
			fgets(line, sizeof(line), fp);
			if(!isspace(*line) && (*line)) 
			{
				if (!QUASI_UBG) 
				{
					sensorPool[count].updateNeighbors(line);
					++count;
				}
				else 
				{
					//sensorPool[i].updateNeighborsQuasiUBG();
				}	
			}
		}while(!feof(fp));
		fclose(fp);
		if(count != nSensors)
		{
			fprintf_s(stderr,"nSenor! Error in %s %d\n",  __FILE__, __LINE__);
			exit(0);
		}
	}
	// Test if the boundary nodes are connected.
	fprintf_s(stderr,"Testing the boundary nodes to see if the boundaries are connected...\n");
	for (int iType = 0; i <= MAX_BOUNDARY_TYPE; ++ iType)
	{
		for ( i=0; i < nSensors; ++ i)
		{
			Sensor* s = &(sensorPool[i]);
			if (s->isBelongToType(iType) )
			{
				// Lamda expression is very useful.
				auto nbIt = find_if(s->neighbors.begin(),s->neighbors.end(),[&, iType] (Sensor* lamdaS) -> bool
				{
					return ( lamdaS->isBelongToType(iType)  );
				});

				if (s->neighbors.end() == nbIt)
				{
					fprintf_s(stderr,"An boundary node have no boundary neighbors (same boundary.)! Error in %s %d\n",  __FILE__, __LINE__);
					exit(0);
				}
			}
		}
	}
	fprintf_s(stderr,"Boundary sensors are connected.\n");

	fprintf_s(stderr,"Calculating the degree...\n");

	double totalDist = 0;
	int totalHops = 0;		
	double totalDeg = 0;
	double totalDeg2 = 0;

	double totalOverallNodeDegree = 0;

	vector<int> degrees;
	vector<int> degreesN;
	int N = 2;
	int idx1 = 0, idx2 = 0;


	for(/*int*/ i = 0; i < nSensors; ++i) 
	{
		Sensor &s = sensorPool[i];
		totalOverallNodeDegree += s.neighbors.size();
	}


	UniformRandVar  samplerand;
	//	int i;
	for(i = 0; i < nSensors; ++i) 
	{
		if (samplerand.value() < 0.1) 
		{
			degrees.push_back(sensorPool[i].neighbors.size());
			degreesN.push_back(sensorPool[i].getMultiHopNeighborhoodSize(N));

		}
	}

	make_heap(degrees.begin(), degrees.end());
	idx1 = int(0.2 * degrees.size());
	idx2 = int(0.3 * degrees.size());

	for(i = idx1; i < idx2; ++i) 
	{
		totalDeg += degrees[i];
	}

	double avg1 = totalDeg / double(idx2 - idx1);

	make_heap(degreesN.begin(), degreesN.end());
	idx1 = int(0.2 * degreesN.size());
	idx2 = int(0.2*N * degreesN.size());

	for(i = idx1; i < idx2; ++i) 
	{
		totalDeg2 += degreesN[i];
	}

	double avgN = totalDeg2 / double(idx2 - idx1);

	fprintf_s(stderr,"measured: avg1hop = %g, %g\n", avg1, avgN);

	avg1HopNeighborhoodSize = avg1;
	avgNHopNeighborhoodSize = avgN;
	avgOverallNodeDegree = totalOverallNodeDegree/nSensors;
	//High_Edge_Threshold = T_NS/avgNHopNeighborhoodSize;
}
void Field::calSensorImportance()
{
	int sNHopNeighborSize = 0;
	const int NHOP_SCALE = 2;
	for (int i = 0; i < nSensors; ++ i)
	{
		Sensor* s = &sensorPool[i];
		if (s->boundarySensorType >= 0)
		{
			// Calculate n hop neighborhood size.
			sNHopNeighborSize = s->getMultiHopNeighborhoodSize(NHOP_SCALE);
			// Compute criticality.
			s->Criticality = (double)sNHopNeighborSize/(avgNHopNeighborhoodSize/2) - 1.0;
			// Calculate node importance.
			s->Rou_p = fabs(s->Criticality) + 1;
		}
		
	}

}

Sensor* Field::findMostImportantNode(const int iType, vector<Sensor*> & boundarySensors)
{
	
	double highestImportance = -1;// importance is a positive number.
	Sensor* highestImportanceSensor = NULL;
	//assert(!boundarySensors.empty());
	if (boundarySensors.empty())
	{
		return highestImportanceSensor;
	}
	//else the flooding did not cover all the boundary nodes.
	for (auto it = boundarySensors.begin(); it != boundarySensors.end(); ++ it)
	{
		Sensor* s = *it;
		if (s->iType2inactive[iType] || s->iType2Landmark[iType])// a sensor is chosen as a landmark or is activated by a landmark.
		{
			continue;
		}
		else
		{
			if (s->Rou_p > highestImportance)
			{
				highestImportance = s->Rou_p;
				highestImportanceSensor = s;
			}
		}
	}
	return highestImportanceSensor;
}
void Field::determineEdgeNodes()
{
	int confirm_num = 0;
	int possible_num = 0;
	for(int i = 0; i < nSensors; ++i) 
	{
		Sensor* s = sensorPool + i;
		
		int SNNeighborhoodsize = s->getMultiHopNeighborhoodSize(2);
		
		//Compute criticality
		s->Criticality = (double)SNNeighborhoodsize/(avgNHopNeighborhoodSize/2) - 1.0;
		if(s->Criticality <= Convex_threshold)
		{
			s->convex = true;
		}
		if ((double)SNNeighborhoodsize <= Low_Edge_Threshold * avgNHopNeighborhoodSize) 
		{
			s->onEdge_Degree = 1.0;							//the smaller, the closer to the middle
			++confirm_num;
		}
		else if((double)SNNeighborhoodsize > Low_Edge_Threshold * avgNHopNeighborhoodSize && (double)SNNeighborhoodsize < High_Edge_Threshold * avgNHopNeighborhoodSize)
		{
			double k = 1.0/(High_Edge_Threshold - Low_Edge_Threshold);
			s->onEdge_Degree = -k * (SNNeighborhoodsize/avgNHopNeighborhoodSize) + High_Edge_Threshold * k;
			++possible_num;
		}
		else
			s->onEdge_Degree = 0;
	}

	int count = 0;
	SCITER neighborpos;

	for(int i = 0; i < nSensors; ++i)
	{
		
		printf("Processing the %dth sensor.\n", i);
		
		if(sensorPool[i].onEdge_Degree > 0)
		{
			BoundarySensorSet.push_back(&sensorPool[i]);
			count++;
		}
		if(count > confirm_num + possible_num)
		{
			fprintf_s(stderr, "Number of edge points incorrect!");
			exit(0);
		}
	}
	printf("The confirmed number of edge: %d, The possible number of edge: %d\n", confirm_num, possible_num);
}

void Field::refineEdgeNodes()
{
	unsigned int i;
	int within_edge_points = 0;
	int edgenodeset_size = BoundarySensorSet.size();
	for(i = 0; i < BoundarySensorSet.size(); ++i)
	{
		Sensor *b = BoundarySensorSet[i];
		within_edge_points = b->LocalFloodingForRefineBoundaryPoints(2);
		//if(within_edge_points <= avg1HopNeighborhoodSize/3*2)
		if(within_edge_points <= avg1HopNeighborhoodSize/3)
		{
			b->onEdge_Degree = 0;
			printf("Cancel this boundary sensor %d\n", b->index);
		}
	}

	BoundarySensorSet.clear();
	for(i = 0; i < nSensors; ++i)
	{
		if(sensorPool[i].onEdge_Degree != 0)
		{
			BoundarySensorSet.push_back(&sensorPool[i]);
		}
	}
	printf("Initial number of boundary points %d, Total number of boundary points %d\n", edgenodeset_size, BoundarySensorSet.size());
}

void Field::refineConvexpoints(int khop, int within_num)
{
	vector<Sensor*> convexset;
	int convexset_size = 0;
	int within_convex_points = 0;
	for(unsigned int i = 0; i < BoundarySensorSet.size(); ++i)
	{
		if(BoundarySensorSet[i]->convex)
		{
			convexset.push_back(BoundarySensorSet[i]);
		}
	}
	convexset_size = convexset.size();
	printf("Total number of convex points %d\n", convexset_size);
	for(unsigned int i = 0; i < convexset.size(); ++i)
	{
		Sensor *c = convexset[i];
		within_convex_points = c->LocalFloodingForRefineConvexPoints(khop);
		if(within_convex_points <= within_num)
		{
			c->convex = false;
			--convexset_size;
			printf("Cancel this convex sensor %d\n", c->index);
		}
	}
	printf("Total number of convex points %d\n", convexset_size);
}

void Field::refineConcavepoints(int khop, int within_num)
{
	vector<Sensor*> concaveset;
	int concaveset_size = 0;
	int within_concave_points = 0;
	for(unsigned int i = 0; i < BoundarySensorSet.size(); ++i)
	{
		if(BoundarySensorSet[i]->concave)
		{
			concaveset.push_back(BoundarySensorSet[i]);
		}
	}
	concaveset_size = concaveset.size();
	printf("Total number of concave points %d\n", concaveset_size);
	for(unsigned int i = 0; i < concaveset.size(); ++i)
	{
		Sensor *c = concaveset[i];
		within_concave_points = c->LocalFloodingForRefineConcavePoints(khop);
		if(within_concave_points <= within_num)
		{
			c->concave = false;
			--concaveset_size;
			printf("Cancel this concave sensor %d\n", c->index);
		}
	}
	printf("Total number of concave points %d\n", concaveset_size);
}
//for chicago airport
#if 0
void Field::detectConcavePoint(int rhop, int dealingpoints)
{
	unsigned int i;
	double lowestCRI = 0xffff;
	double biggestCRI = 0;
	double AvgCriticality_Threshold = 0;
	Sensor *CenterSensor = NULL;
	int number_of_withincircle = 0;
	int localcovered_number = 0;
	int current_points = 0;
	//get the biggest criticality point to be the start point

	for( i = 0; i < BoundarySensorSet.size(); ++i)
	{
		if(biggestCRI < BoundarySensorSet[i]->Criticality)
		{
			biggestCRI = BoundarySensorSet[i]->Criticality;
			CenterSensor = BoundarySensorSet[i];
		}
	}

	//get the random point to be the start point
	//UniformRandVar urandv(0, BoundarySensorSet.size()-0.1);
	//int random_value = (int)urandv.value();
	//CenterSensor = BoundarySensorSet[random_value];
	Sensor *PathStart = NULL;
	Sensor *PathEnd = NULL;
	Path DetectConcavePointPath;

	//start to establish circle areas
	vector<CircleArea> OriginalAreaSet;
	vector<CircleArea> FinalAreaSet;
	for(unsigned int j = 0; j < BoundarySensorSet.size(); ++j)
	{
		if(BoundarySensorSet[j]->within_circle)
		{
			++number_of_withincircle;
		}
	}
	printf("From the beginning, number of withincircle is %d\n", number_of_withincircle);
	while (number_of_withincircle != BoundarySensorSet.size())
	{
		number_of_withincircle = 0;


		CircleArea CA(this, CenterSensor, rhop);

		printf("There are %d convex points contain in the circle area!\n", CA.contain_convex_number);
		//reselect center
		//lowestCRI = 0xffff;
		biggestCRI = -100;
		for(unsigned int j = 0; j < BoundarySensorSet.size(); ++j)
		{
			if(BoundarySensorSet[j]->within_circle)
			{
				++number_of_withincircle;
			}
			else
			{
				if(biggestCRI < BoundarySensorSet[j]->Criticality)
				{
					biggestCRI = BoundarySensorSet[j]->Criticality;
					CenterSensor = BoundarySensorSet[j];
				}
			}
		}
		/*
		random_value = (int)urandv.value();
		while(BoundarySensorSet[random_value]->within_circle)
		{
		random_value = (int)urandv.value();
		}
		CenterSensor = BoundarySensorSet[random_value];
		*/
		if(!CA.Circle.empty() && CA.contain_convex_number <= Contain_convex_threshold /*&& CA.AreaSensorSet.size() > avgNHopNeighborhoodSize*/)
		{
			OriginalAreaSet.push_back(CA);
			//dumpCircle("Circle", CA);
		}
		for(int j = 0; j < nSensors; ++j)
		{
			sensorPool[j].broadcastCovered = 0;
			sensorPool[j].on_boundary_circle = false;
		}
		printf("%d sensors are within circle!\n", number_of_withincircle);
	}

	//Picking efficient circles
	double sum_avgCRI = 0;
	for(i = 0; i < OriginalAreaSet.size(); ++i)
	{
		sum_avgCRI += OriginalAreaSet[i].AvgCriticality;
	}
	AvgCriticality_Threshold = sum_avgCRI/OriginalAreaSet.size();
	for (i = 0; i < OriginalAreaSet.size(); ++i)
	{
		if(OriginalAreaSet[i].AvgCriticality >= RATIO*AvgCriticality_Threshold)
		{
			FinalAreaSet.push_back(OriginalAreaSet[i]);
			//dumpCircle("Circle", OriginalAreaSet[i]);
		}
	}

	//deal every circles
	int concave_points = 0;
	int dump_path = 0;
	set<Sensor*> DealingSensors;
	for(i = 0; i < FinalAreaSet.size(); ++i)
	{
		DealingSensors.clear();
		dump_path = 0;
		printf("Processing area %d, total areas %d!\n", i, FinalAreaSet.size());
		//if(FinalAreaSet[i].Center->index == 19777)
		//	printf("wait!\n");
		current_points = 0;
		concave_points = 0;
		multimap<double, Sensor*, less<double> > Circle = FinalAreaSet[i].GetCircle();
		for(multimap<double, Sensor*, less<double> >::iterator CCpos = Circle.begin(); CCpos != Circle.end(); ++CCpos)
		{
			Sensor *d = CCpos->second;
			if(d->Criticality < 0.3)
			{
				DealingSensors.insert(d);
				++current_points;
				if(current_points >= dealingpoints)
				{
					break;
				}
			}

		}
		printf("Need to dealing %d sensors\n", DealingSensors.size());
		DetectConcavePointPath.pafield = this;
		set<Sensor*> concave_in_circle;
		concave_in_circle.clear();
		for(SCITER DSpos = DealingSensors.begin(); DSpos != DealingSensors.end(); ++DSpos)
		{
			PathStart = *DSpos;

			for(SCITER newDspos = DSpos; newDspos != DealingSensors.end(); ++newDspos)
			{
				PathEnd = *newDspos;
				if(PathEnd != PathStart)
				{
					//DetectConcavePointPath.EstablishPathBySideStep(PathStart, PathEnd, FinalAreaSet[i].AreaSensorSet);
					DetectConcavePointPath.EstablishLowestCriticalityShortestPath(PathStart, PathEnd, FinalAreaSet[i].AreaSensorSet);

					for(unsigned int k = 0; k < DetectConcavePointPath.path.size(); ++k)
					{
						Sensor *dp = DetectConcavePointPath.path[k];
						if(dp->Criticality <= Concave_threshold && dp->onEdge_Degree != 0 && dp->WithinArea(FinalAreaSet[i].AreaSensorSet))
						{
							dp->concave = true;
							concave_in_circle.insert(dp);
							++dump_path;
						}
					}
					if( i == 1 )
					{
						dumpPath("Path", DetectConcavePointPath);
					}
				}
			}
		}
		if(concave_in_circle.size() < 0)//rhop)
		{
			printf("This circle may not on concave area, cancel all the concave point in it!\n");
			for(SCITER FASpos = FinalAreaSet[i].AreaSensorSet.begin(); FASpos != FinalAreaSet[i].AreaSensorSet.end(); ++FASpos)
			{
				Sensor *f = *FASpos;
				if(f->concave)
				{
					f->concave = false;
				}
			}
		}
		for(unsigned int m = 0; m < BoundarySensorSet.size(); ++m)
		{
			if(BoundarySensorSet[m]->concave)
			{
				++concave_points;
			}
		}
		printf("Detect %d concave points!\n", concave_points);

	}
}
#endif
//for cube with H
void Field::detectConcavePoint(int rhop, int dealingpoints)
{
	unsigned int i;
	double lowestCRI = 0xffff;
	double biggestCRI = 0;
	double AvgCriticality_Threshold = 0;
	Sensor *CenterSensor = NULL;
	int number_of_withincircle = 0;
	int localcovered_number = 0;
	int current_points = 0;
	//get the biggest criticality point to be the start point
	
	for( i = 0; i < BoundarySensorSet.size(); ++i)
	{
		if(biggestCRI < BoundarySensorSet[i]->Criticality)
		{
			biggestCRI = BoundarySensorSet[i]->Criticality;
			CenterSensor = BoundarySensorSet[i];
		}
	}
	
	//get the random point to be the start point
	//UniformRandVar urandv(0, BoundarySensorSet.size()-0.1);
	//int random_value = (int)urandv.value();
	//CenterSensor = BoundarySensorSet[random_value];
	Sensor *PathStart = NULL;
	Sensor *PathEnd = NULL;
	Path DetectConcavePointPath;

	//start to establish circle areas
	vector<CircleArea> OriginalAreaSet;
	vector<CircleArea> FinalAreaSet;
	for(unsigned int j = 0; j < BoundarySensorSet.size(); ++j)
	{
		if(BoundarySensorSet[j]->within_circle)
		{
			++number_of_withincircle;
		}
	}
	printf("From the beginning, number of withincircle is %d\n", number_of_withincircle);
	while (number_of_withincircle != BoundarySensorSet.size())
	{
		number_of_withincircle = 0;
		
		
		CircleArea CA(this, CenterSensor, rhop);
		
		printf("There are %d convex points contain in the circle area!\n", CA.contain_convex_number);
		//reselect center
		//lowestCRI = 0xffff;
		biggestCRI = -100;
		for(unsigned int j = 0; j < BoundarySensorSet.size(); ++j)
		{
			if(BoundarySensorSet[j]->within_circle)
			{
				++number_of_withincircle;
			}
			else
			{
				if(biggestCRI < BoundarySensorSet[j]->Criticality)
				{
					biggestCRI = BoundarySensorSet[j]->Criticality;
					CenterSensor = BoundarySensorSet[j];
				}
			}
		}
		/*
		random_value = (int)urandv.value();
		while(BoundarySensorSet[random_value]->within_circle)
		{
			random_value = (int)urandv.value();
		}
		CenterSensor = BoundarySensorSet[random_value];
		*/
		if(!CA.Circle.empty() && CA.contain_convex_number <= Contain_convex_threshold/* && CA.AreaSensorSet.size() > avgNHopNeighborhoodSize*/)
		{
			OriginalAreaSet.push_back(CA);
			//dumpCircle("Circle", CA);
		}
		for(int j = 0; j < nSensors; ++j)
		{
			sensorPool[j].broadcastCovered = 0;
			sensorPool[j].on_boundary_circle = false;
		}
		printf("%d sensors are within circle!\n", number_of_withincircle);
	}
	
	//Picking efficient circles
	double sum_avgCRI = 0;
	for(i = 0; i < OriginalAreaSet.size(); ++i)
	{
		sum_avgCRI += OriginalAreaSet[i].AvgCriticality;
	}
	AvgCriticality_Threshold = sum_avgCRI/OriginalAreaSet.size();
	for (i = 0; i < OriginalAreaSet.size(); ++i)
	{
		if(OriginalAreaSet[i].AvgCriticality >= RATIO*AvgCriticality_Threshold)
		{
			FinalAreaSet.push_back(OriginalAreaSet[i]);
			//dumpCircle("Circle", OriginalAreaSet[i]);
		}
	}
	
	//deal every circles
	int concave_points = 0;
	int dump_path = 0;
	set<Sensor*> DealingSensors;
	for(i = 0; i < FinalAreaSet.size(); ++i)
	{
		DealingSensors.clear();
		dump_path = 0;
		printf("Processing area %d, total areas %d!\n", i, FinalAreaSet.size());
		//if(FinalAreaSet[i].Center->index == 19777)
		//	printf("wait!\n");
		current_points = 0;
		concave_points = 0;
		multimap<double, Sensor*, less<double> > Circle = FinalAreaSet[i].GetCircle();
		for(multimap<double, Sensor*, less<double> >::iterator CCpos = Circle.begin(); CCpos != Circle.end(); ++CCpos)
		{
			Sensor *d = CCpos->second;
			//if(d->Criticality < 0.5)
			//{
				DealingSensors.insert(d);
				++current_points;
				if(current_points >= Circle.size()/4*3)
				{
					break;
				}
			//}
		}
		printf("Need to dealing %d sensors\n", DealingSensors.size());
		DetectConcavePointPath.pafield = this;
		set<Sensor*> concave_in_circle;
		concave_in_circle.clear();
		for(SCITER DSpos = DealingSensors.begin(); DSpos != DealingSensors.end(); ++DSpos)
		{
			PathStart = *DSpos;
			
			for(SCITER newDspos = DSpos; newDspos != DealingSensors.end(); ++newDspos)
			{
				PathEnd = *newDspos;
				if(PathEnd != PathStart)
				{
					//DetectConcavePointPath.EstablishPathBySideStep(PathStart, PathEnd, FinalAreaSet[i].AreaSensorSet);
					DetectConcavePointPath.EstablishLowestCriticalityShortestPath(PathStart, PathEnd, FinalAreaSet[i].AreaSensorSet);
					
					for(unsigned int k = 0; k < DetectConcavePointPath.path.size(); ++k)
					{
						Sensor *dp = DetectConcavePointPath.path[k];
						if(dp->Criticality > Concave_threshold && dp->onEdge_Degree != 0 && dp->WithinArea(FinalAreaSet[i].AreaSensorSet))
						{
							dp->concave = true;
							concave_in_circle.insert(dp);
							++dump_path;
						}
					}
					if( i == 1 )
					{
						dumpPath("Path", DetectConcavePointPath);
					}
				}
			}
		}
		if(concave_in_circle.size() < 0)//rhop)
		{
			printf("This circle may not on concave area, cancel all the concave point in it!\n");
			for(SCITER FASpos = FinalAreaSet[i].AreaSensorSet.begin(); FASpos != FinalAreaSet[i].AreaSensorSet.end(); ++FASpos)
			{
				Sensor *f = *FASpos;
				if(f->concave)
				{
					f->concave = false;
				}
			}
		}
		for(unsigned int m = 0; m < BoundarySensorSet.size(); ++m)
		{
			if(BoundarySensorSet[m]->concave)
			{
				++concave_points;
			}
		}
		printf("Detect %d concave points!\n", concave_points);
		
	}
}

void Field::detectSaddlepoints(int detect_dist, int radius)
{
	ComputeFieldAvghopDistWithinCircle(radius);
	//avgOverallHopDistWithinCircle = 6.546867;
	int i;
	int BoundarySensorSet_size = BoundarySensorSet.size();
	set<Sensor*> convex_set;
	set<Sensor*> concave_set;
	set<Sensor*> convex_boundary_set;//the boundary nodes within convex set
	for(i = 0; i < BoundarySensorSet_size; ++i)
	{
		BoundarySensorSet[i]->asCenter = false;
		if(BoundarySensorSet[i]->convex)
		{
			convex_set.insert(BoundarySensorSet[i]);
		}
	}
	//Picking convex node in convex boundary
	//Get avg convex neighbor size
	for(SCITER cspos1 = convex_set.begin(); cspos1 != convex_set.end(); ++cspos1)
	{
		Sensor *s1 = *cspos1;
		for(SCITER cspos2 = convex_set.begin(); cspos2 != convex_set.end(); ++cspos2)
		{
			Sensor *s2 = *cspos2;
			if(s1 != s2 && s1->location.distance(s2->location) <= COMM_RANGE)
			{
				s1->convex_neighbor_size++;
			}
		}
	}
	int sum_convex_neighbor_size = 0;
	double avg_convex_neighbor_size = 0;
	for(SCITER cspos = convex_set.begin(); cspos != convex_set.end(); ++cspos)
	{
		Sensor *s = *cspos;
		sum_convex_neighbor_size += s->convex_neighbor_size;
	}
	avg_convex_neighbor_size = (double)sum_convex_neighbor_size/convex_set.size();
	//Picking
	for(SCITER cspos = convex_set.begin(); cspos != convex_set.end(); ++cspos)
	{
		Sensor *s = *cspos;
		if(s->convex_neighbor_size < 2)
		{
			convex_boundary_set.insert(s);
		}
	}
	SCITER cvs_pos;
	int saddle_number = 0;
	for(cvs_pos = convex_boundary_set.begin(); cvs_pos != convex_boundary_set.end(); ++cvs_pos)
	{
		Sensor *cv = *cvs_pos;
		if(cv->FloodingtoFindConcavepointsWithin(detect_dist, concave_set))
		{
			/*
			for(SCITER ccs_pos = concave_set.begin(); ccs_pos != concave_set.end(); ++ccs_pos)
			{
				Sensor *cc = *ccs_pos;
				Path cp(cv, cc);
				cp.EstablishPath(cv, cc);
				for(unsigned int j = 1; j < cp.path.size()-1; ++j)
				{
					Sensor *center_circle = cp.path[j];
					if(!center_circle->asCenter)
					{
						CircleArea C_for_Saddle(this, center_circle, radius);
						C_for_Saddle.ComputingCircleAvgHopDist();
						printf("Circle center %d->avgHopDist %lf, Field overall avghopdist %lf\n", center_circle->index, C_for_Saddle.AvgHopdist, avgOverallHopDistWithinCircle);
						if( (C_for_Saddle.AvgHopdist-avgOverallHopDistWithinCircle)/avgOverallHopDistWithinCircle >= Saddle_threshold )
						{
							C_for_Saddle.Center->saddle = true;
							++saddle_number;
							printf("Detect a saddle point!\n");
						}
						for(int k = 0; k < BoundarySensorSet_size; ++k)
						{
							BoundarySensorSet[k]->within_circle = false;
						}
					}
				}
			}
			*/
			set<Sensor*> DealingSet;
			set<Sensor*> FinalDealingSet;
			/*
			for(SCITER cvneighborspos = cv->neighbors.begin(); cvneighborspos != cv->neighbors.end(); ++cvneighborspos)
			{
				Sensor *v = *cvneighborspos;
				if(!v->convex && !v->concave && v->onEdge_Degree != 0)
				{
					DealingSet.insert(v);
				}
			}
			*/
			cv->getMultiHopNeighborhood(2, DealingSet);
			int DealingNum = 0;
			double sum_CN = 0;
			double avg_CN = 0;
			for(SCITER DealingSetpos = DealingSet.begin(); DealingSetpos != DealingSet.end(); ++DealingSetpos)
			{
				Sensor *de = *DealingSetpos;
				if(!de->convex && !de->concave && de->onEdge_Degree != 0)
				{
					de->KHOPConvexNeighborSize(2);
					++DealingNum;
					sum_CN += de->khopConvexNeighborSize;
				}
			}
			avg_CN = sum_CN/(double)DealingNum;
			for(SCITER DealingSetpos = DealingSet.begin(); DealingSetpos != DealingSet.end(); ++DealingSetpos)
			{
				Sensor *de = *DealingSetpos;
				if(!de->convex && !de->concave && de->onEdge_Degree != 0 && de->khopConvexNeighborSize < avg_CN)
				{
					if(de->Criticality >= Concave_threshold1 && de->Criticality <= Concave_threshold2)
					{
						de->concave = true;
					}
					else
					{
						FinalDealingSet.insert(de);
					}
				}
			}
			for(int k = 0; k < BoundarySensorSet_size; ++k)
			{
				BoundarySensorSet[k]->within_circle = false;
			}
			for(SCITER Dspos = FinalDealingSet.begin(); Dspos != FinalDealingSet.end(); ++Dspos)
			{
				Sensor *d = *Dspos;
				CircleArea C_for_Saddle(this, d, radius);
				C_for_Saddle.ComputingCircleAvgHopDist();
				dumpCircle("Circle", C_for_Saddle);
				printf("Circle center %d->avgHopDist %lf, Field overall avghopdist %lf\n", d->index, C_for_Saddle.AvgHopdist, avgOverallHopDistWithinCircle);
				if( (C_for_Saddle.AvgHopdist-avgOverallHopDistWithinCircle)/avgOverallHopDistWithinCircle >= Saddle_threshold )
				{
					C_for_Saddle.Center->saddle = true;
					++saddle_number;
					printf("Detect a saddle point!\n");
				}
				for(int k = 0; k < BoundarySensorSet_size; ++k)
				{
					BoundarySensorSet[k]->asCenter = false;
					BoundarySensorSet[k]->on_boundary_circle = false;
					BoundarySensorSet[k]->within_circle = false;
				}
			}
		}
	}
	printf("Field overall avg_hop_dist within a circle: %lf\n", avgOverallHopDistWithinCircle);
	printf("There are %d saddle points!\n", saddle_number);
}

void Field::ComputeFieldAvghopDistWithinCircle(int rhop)
{
	unsigned int i;
	//double lowestCRI = 0xffff;
	vector<double> hoplist;
	double sum_hop_dist = 0;
	Sensor *CenterSensor = NULL;
	int sample_times = 30;
	UniformRandVar randv(0, BoundarySensorSet.size());
	//generate circle
	for( i = 0; i < BoundarySensorSet.size(); ++i)
	{
		BoundarySensorSet[i]->within_circle = false;
		BoundarySensorSet[i]->asCenter = false;
		
	}
	int Bindex = 0;
	do
	{
		Bindex = (int)randv.value();
		CenterSensor = BoundarySensorSet[Bindex];
	} while(BoundarySensorSet[Bindex]->convex || BoundarySensorSet[Bindex]->concave);
	int count = 0;
	while (count != sample_times)
	{
		CircleArea CA(this, CenterSensor, rhop);
		//biggestCRI = -100;
		do
		{
			Bindex = (int)randv.value();
			CenterSensor = BoundarySensorSet[Bindex];
		} while(BoundarySensorSet[Bindex]->convex || BoundarySensorSet[Bindex]->concave || BoundarySensorSet[Bindex]->within_circle);
		if(!CA.Circle.empty())
		{
			++count;
			CA.ComputingCircleAvgHopDist();
			hoplist.push_back(CA.AvgHopdist);
			//sum_hop_dist += CA.AvgHopdist;
			printf("Circle center %d, AvgHopDist %lf\n", CA.Center->index, CA.AvgHopdist);
		}
		for(multimap<double, Sensor*, less<double> >::iterator Circlepos = CA.Circle.begin(); Circlepos != CA.Circle.end(); ++Circlepos)
		{
			Sensor *c = Circlepos->second;
			c->KnownHopDistEndID.clear();
			c->ID_HopDist_List.clear();
		}	
	}
	sort(hoplist.begin(), hoplist.end());
	for(unsigned int h = hoplist.size()/5; h < hoplist.size()/5*4; ++h)
	{
		sum_hop_dist += hoplist[h];
	}
	avgOverallHopDistWithinCircle = sum_hop_dist/((double)sample_times/5*3);
	printf("Average overall hop distance within circle: %lf\n", avgOverallHopDistWithinCircle);
}
void Field::selectLandMarks(const int kHOP)
{
	////Hop distance between landmarks is larger than 1.
	//assert(kHOP > 1);
	//vector<Sensor*> tmpBoundarySensors;
	//Sensor* importantNode = NULL;
	//for (int iType=1; iType <= MAX_BOUNDARY_TYPE; ++ iType)
	//{
	//	// Find sensors on each type of boundary.
	//	for (int iSensor = 0; iSensor < nSensors; ++ iSensor )
	//	{
	//		Sensor* s = &sensorPool[iSensor];
	//		if (s->isBelongToType(iType))
	//		{
	//			// At this step, all the boundary sensors are not landmarks and should be not active.
	//			s->landmark = false;
	//			s->inactive = false;
	//			s->broadcastCovered = 0;
	//			tmpBoundarySensors.push_back(s);
	//		}
	//	}
	//	// Initialize important node for the while loop;
	//	importantNode = findMostImportantNode(tmpBoundarySensors);
	//	// Flooding has a boundary.
	//	vector<Sensor*> localFloodingBoundary;
	//	// If can not find a important node on the local flooding boundary, then the landmarks are found.
	//	while(importantNode != NULL)
	//	{
	//		importantNode->landmark = true;
	//		LandmarkSet.push_back(importantNode);
	//		importantNode -> floodingOnBoundary(iType,kHOP,localFloodingBoundary);
	//		importantNode = findMostImportantNode(localFloodingBoundary);			
	//	}
	//	cerr<<"Now there are "<<LandmarkSet.size()<<" landmarks."<<endl;
	//	for_each(tmpBoundarySensors.begin(),tmpBoundarySensors.end(),[&](Sensor* lamdaS)
	//	{
	//		lamdaS->broadcastCovered = 0;	
	//	});

	//	tmpBoundarySensors.clear();
	//}
}
void Field::PickingLandmarkWithSignificance(int khop)
{
	assert(khop > 1);
	vector<Sensor*> tmpBoundarySensors;
	/*LandmarkSet.clear();*/

	for (int iType = 1; iType <= MAX_BOUNDARY_TYPE; ++ iType)
	{
		iType2LandmarkSet[iType].clear();
		// Find sensors on each type of boundary.
		for (int iSensor = 0; iSensor < nSensors; ++ iSensor )
		{
			Sensor* s = &sensorPool[iSensor];
			if (s->isBelongToType(iType))
			{
				// At this step, all the boundary sensors are not landmarks and should be not active.
				s->iType2Landmark[iType] = false;
				s->iType2inactive[iType] = false;
				s->broadcastCovered = 0;
				tmpBoundarySensors.push_back(s);
			}
		}

		Sensor *floodingSensor = NULL;
		int landmark_space = 0;

		//Select a good start sensor whose Rou_p is largest.
		floodingSensor = findMostImportantNode(iType,tmpBoundarySensors);

		int count = 0;
		int non_khop_landmarkspace = 0;


		while(floodingSensor != NULL)
		{
			//double Rp = floodingSensor->Rou_p;

			landmark_space = khop;//(int)((double)khop/Rp+0.5);
			
			floodingSensor->floodingOnBoundary(iType,landmark_space);
			// floodingSensor->landmark = true;
			// floodingSensor->broadcastCovered = 1;
			iType2LandmarkSet[iType].push_back(floodingSensor);
			iType2LandmarkSet[iType].erase(remove_if(iType2LandmarkSet[iType].begin(),iType2LandmarkSet[iType].end(),[&](Sensor* s)
			{
				return (!s->iType2Landmark[iType]);
			}),	iType2LandmarkSet[iType].end());
			floodingSensor = findMostImportantNode(iType,tmpBoundarySensors);
		}

		//////////////////////////////////////////////////////////////////////////
		//Refine landmarks selecting
		bool problemtag;
		do 
		{
			problemtag = false;
			for(vector<Sensor>::size_type i = 0; i < iType2LandmarkSet[iType].size(); ++i)
			{
				/*double Rp = LandmarkSet[i]->Rou_p;*/

				landmark_space = khop;/*(int)((double)khop/Rp+0.5);*/
				int hop = iType2LandmarkSet[iType][i]->distanceBetweenLandmarks(iType,landmark_space);
				// low criticality must be removed?
				if(hop < landmark_space /*|| fabs(LandmarkSet[i]->Criticality-1.0) < 0.1 */)
				{
					iType2LandmarkSet[iType][i]->landmark = false;
					iType2LandmarkSet[iType][i]->inactive = true;
					problemtag = true;
					printf("Landmark selecting has problems, wrong landmark %d, hop %d space %d\n", iType2LandmarkSet[iType][i]->index, hop,landmark_space);
					
					//printf("Total landmarks %d\n", Number_of_Landmark);
					break;
				}
			}
			iType2LandmarkSet[iType].erase(remove_if(iType2LandmarkSet[iType].begin(),iType2LandmarkSet[iType].end(),[&](Sensor* s)
			{
				return (!s->iType2Landmark[iType]);
			}),iType2LandmarkSet[iType].end());

			// printf("Total landmarks %d\n", LandmarkSet.size());
		} while (problemtag);
		//detecting landmark state
		for(vector<Sensor>::size_type i = 0; i < tmpBoundarySensors.size(); ++i)
		{
			if( (    tmpBoundarySensors[i]->iType2Landmark[iType] 
				&&   tmpBoundarySensors[i]->iType2inactive[iType]) 
				|| (!tmpBoundarySensors[i]->iType2Landmark[iType]  && !tmpBoundarySensors[i]->iType2inactive[iType] ) )
			{
				fprintf_s(stderr,"landmark states are both landmark and inactive!!!!!\n");
				exit(EXIT_FAILURE);
			}
		}
	
		cerr<<"The "<<iType<<" level has "<<iType2LandmarkSet[iType].size()<<" landmarks."<<endl;

		for_each(tmpBoundarySensors.begin(),tmpBoundarySensors.end(),[&](Sensor* lamdaS)
		{
			lamdaS->broadcastCovered = 0;	
		});

		tmpBoundarySensors.clear();
	}
	
}

#if 0
void Field::SelectingLamdmarks(int khop)
{
	Sensor *floodingSensor = NULL;
	unsigned int i = 0;

	//Select a good start sensor whose OnEdgeDegree is smallest
	floodingSensor = EdgeSensorSet[0];
	for(i = 1; i < EdgeSensorSet.size(); i++)
	{
		if(floodingSensor->onEdge_Degree < EdgeSensorSet[i]->onEdge_Degree)
			floodingSensor = EdgeSensorSet[i];
	}
	floodingSensor->broadcastCovered = 1;
	floodingSensor->landmark = true;
	++Number_of_Landmark;
	int count = 0;
	do 
	{
		if(floodingSensor == NULL)
		{
			floodingSensor = ReselectingFloodingSensor();
			count = 0;
		}
		while(floodingSensor != NULL)
		{
			floodingSensor = floodingSensor->LocalFlooding(khop);
		}
		printf("Total landmarks %d\n", Number_of_Landmark);

		for(i = 0; i < EdgeSensorSet.size(); ++i)
		{
			if(EdgeSensorSet[i]->inactive || EdgeSensorSet[i]->landmark)
				++count;
		}
		printf("Covered sensors %d, Edge sensors %d\n",count, EdgeSensorSet.size());
	} while (count != EdgeSensorSet.size());
	for(i = 0; i < EdgeSensorSet.size(); ++i)
	{
		EdgeSensorSet[i]->broadcastCovered = 0;
		if(EdgeSensorSet[i]->landmark)
		{
			LandmarkSet.push_back(EdgeSensorSet[i]);
		}
	}
}
#endif

void Field::SelectingLamdmarksRandomly(int khop)
{
	Sensor *floodingSensor = NULL;
	unsigned int i = 0;

	//Select a good start sensor whose OnEdgeDegree is biggest
	floodingSensor = BoundarySensorSet[0];
	for(i = 1; i < BoundarySensorSet.size(); i++)
	{
		if(floodingSensor->onEdge_Degree < BoundarySensorSet[i]->onEdge_Degree)
			floodingSensor = BoundarySensorSet[i];
	}
	floodingSensor->broadcastCovered = 1;
	floodingSensor->landmark = true;
	++Number_of_Landmark;
	int count = 0;
	do 
	{
		while(floodingSensor != NULL)
		{
			floodingSensor->LocalFlooding(khop);
			floodingSensor = ReselectingFloodingSensor();
		}
		printf("Total landmarks %d\n", Number_of_Landmark);
		count = 0;
		for(i = 0; i < BoundarySensorSet.size(); ++i)
		{
			if(BoundarySensorSet[i]->inactive || BoundarySensorSet[i]->landmark)
				++count;
		}
		printf("Covered sensors %d, Edge sensors %d\n",count, BoundarySensorSet.size());
	} while (count != BoundarySensorSet.size());
	for(i = 0; i < BoundarySensorSet.size(); ++i)
	{
		BoundarySensorSet[i]->broadcastCovered = 0;
		/*
		if(EdgeSensorSet[i]->onEdge_Degree < 0.1 && EdgeSensorSet[i]->landmark)
		{
			EdgeSensorSet[i]->landmark = false;
			EdgeSensorSet[i]->inactive = true;
			--Number_of_Landmark;
			printf("onEdge_Degree is too low!!\n");
		}
		*/
		if(BoundarySensorSet[i]->landmark)
		{
			LandmarkSet.push_back(BoundarySensorSet[i]);
		}
	}
	//////////////////////////////////////////////////////////////////////////
	//Refine landmarks selecting
	bool problemtag = false;
	do 
	{
		problemtag = false;
		for(i = 0; i < LandmarkSet.size(); ++i)
		{
			int hop = LandmarkSet[i]->LandmarkFlooding(khop);
			if(hop < khop)
			{
				LandmarkSet[i]->landmark = false;
				LandmarkSet[i]->inactive = true;
				problemtag = true;
				printf("Landmark selecting has problems, wrong landmark %d, hop %d\n", LandmarkSet[i]->index, hop);
				--Number_of_Landmark;
				printf("Total landmarks %d\n", Number_of_Landmark);
				break;
			}
		}
		LandmarkSet.clear();
		for(i = 0; i < BoundarySensorSet.size(); ++i)
		{
			BoundarySensorSet[i]->broadcastCovered = 0;
			if(BoundarySensorSet[i]->landmark)
			{
				LandmarkSet.push_back(BoundarySensorSet[i]);
			}
		}
		printf("Total landmarks %d\n", LandmarkSet.size());
	} while (problemtag);
	//detecting landmark state
	for(i = 0; i < BoundarySensorSet.size(); ++i)
	{
		if( (BoundarySensorSet[i]->landmark && BoundarySensorSet[i]->inactive) || (!BoundarySensorSet[i]->landmark && !BoundarySensorSet[i]->inactive) )
		{
			printf("landmark states are both landmark and inactive!!!!!\n");
			exit(0);
		}
	}
	int countforbadlandmark = 0;
	for( i = 0; i < LandmarkSet.size(); ++i )
	{
		Sensor *s = LandmarkSet[i];
		if(s->onEdge_Degree < 0.1)
		{
			++countforbadlandmark;

		}
	}
	printf("Bad landmark: %d\n", countforbadlandmark);
}

void Field::SelectingLamdmarksRandomlyForSkeleton(int khop)
{
	Number_of_Landmark = 0;
	LandmarkSet.clear();
	EdgeSet.clear();
	TriangleSet.clear();
	HoleSet.clear();
	RefineEdgeSet.clear();
	RefineSensorSet.clear();
	VoronoiCellSensorSet.clear();
	ConnectingLandmarkIDPairs.clear();
	Sensor *floodingSensor = NULL;
	unsigned int i = 0;
	for(i = 0; i < nSensors; ++i)
	{
		sensorPool[i].clear();
		if(sensorPool[i].asSkeleton)
		{
			SkeletonNodeSet.push_back(&sensorPool[i]);
		}
	}
	UniformRandVar randv(0, SkeletonNodeSet.size()-1);
	int rvalue = (int)randv.value();
	//Select a good start sensor whose OnEdgeDegree is biggest
	floodingSensor = SkeletonNodeSet[rvalue];
	
	floodingSensor->broadcastCovered = 1;
	floodingSensor->landmark = true;
	++Number_of_Landmark;
	int count = 0;
	int repeat_times = 0;
	//do 
	//{
	while(floodingSensor != NULL)
	{
		floodingSensor->LocalFloodingforSkeleton(khop);
		int trying = 1000000;
		while(trying != 0)
		{
			rvalue = (int) randv.value();
			floodingSensor = SkeletonNodeSet[rvalue];
			if(!floodingSensor->inactive && !floodingSensor->landmark)
			{
				break;
			}
			--trying;
		}
		if(trying == 0)
		{
			floodingSensor = NULL;
			++repeat_times;
		}
		if(floodingSensor != NULL)
		{
			repeat_times = 0;
			floodingSensor->landmark = true;
			floodingSensor->broadcastCovered = 1;
			++Number_of_Landmark;
			printf("The %dth landmark: sensor %d\n", Number_of_Landmark, floodingSensor->index);
		}
	}
	printf("Total landmarks %d\n", Number_of_Landmark);
	count = 0;
	for(i = 0; i < SkeletonNodeSet.size(); ++i)
	{
		if(SkeletonNodeSet[i]->inactive || SkeletonNodeSet[i]->landmark)
			++count;
	}
	printf("Covered sensors %d, Edge sensors %d\n",count, SkeletonNodeSet.size());
	//} while (count != BoundarySensorSet.size() && repeat_times < 100);
	for(i = 0; i < SkeletonNodeSet.size(); ++i)
	{
		SkeletonNodeSet[i]->broadcastCovered = 0;
		if(SkeletonNodeSet[i]->landmark)
		{
			LandmarkSet.push_back(SkeletonNodeSet[i]);
		}
	}
	//////////////////////////////////////////////////////////////////////////
	//Refine landmarks selecting
	bool problemtag = false;
	do 
	{
		problemtag = false;
		for(i = 0; i < LandmarkSet.size(); ++i)
		{
			int hop = LandmarkSet[i]->LandmarkFlooding(khop);
			if(hop < khop)
			{
				LandmarkSet[i]->landmark = false;
				LandmarkSet[i]->inactive = true;
				problemtag = true;
				printf("Landmark selecting has problems, wrong landmark %d, hop %d\n", LandmarkSet[i]->index, hop);
				--Number_of_Landmark;
				printf("Total landmarks %d\n", Number_of_Landmark);
				break;
			}
		}
		LandmarkSet.clear();
		for(i = 0; i < SkeletonNodeSet.size(); ++i)
		{
			SkeletonNodeSet[i]->broadcastCovered = 0;
			if(SkeletonNodeSet[i]->landmark)
			{
				LandmarkSet.push_back(SkeletonNodeSet[i]);
			}
		}
		printf("Total landmarks %d\n", LandmarkSet.size());
	} while (problemtag);
	//detecting landmark state
	for(i = 0; i < SkeletonNodeSet.size(); ++i)
	{
		if( (SkeletonNodeSet[i]->landmark && SkeletonNodeSet[i]->inactive) )
		{
			printf("landmark states are both landmark and inactive!!!!!\n");
			exit(0);
		}
	}
	int countforbadlandmark = 0;
	for( i = 0; i < LandmarkSet.size(); ++i )
	{
		Sensor *s = LandmarkSet[i];
		if(s->onEdge_Degree < 0.1)
		{
			++countforbadlandmark;

		}
	}
	printf("Bad landmark: %d\n", countforbadlandmark);
}

#if 0
Sensor* Field::ReselectingFloodingSensor()
{
	unsigned int i;
	double maxOEDegree = 0;
	Sensor *result = NULL;
	for(i = 0; i < BoundarySensorSet.size(); ++i)
	{
		if(!BoundarySensorSet[i]->inactive && !BoundarySensorSet[i]->landmark && BoundarySensorSet[i]->onEdge_Degree > 0.1)
		{
			result = BoundarySensorSet[i];
		}
	}
	if(result != NULL)
	{
		result->landmark = true;
		result->broadcastCovered = 1;
		++Number_of_Landmark;
		printf("The %dth landmark: sensor %d\n", Number_of_Landmark, result->index);
	}
	else
	{
		for(i = 0; i < BoundarySensorSet.size(); ++i)
		{
			if(!BoundarySensorSet[i]->inactive && !BoundarySensorSet[i]->landmark)
			{
				result = BoundarySensorSet[i];
			}
		}
		if(result != NULL)
		{
			result->landmark = true;
			result->broadcastCovered = 1;
			++Number_of_Landmark;
			printf("The %dth landmark: sensor %d\n", Number_of_Landmark, result->index);
		}
	}
	
	return result;
}
#endif
Sensor* Field::ReselectingFloodingSensor(vector<Sensor*> & iTypeBoundarySensors)
{
	int itype = 0;
	Sensor *result = NULL;
	result = findMostImportantNode(itype, iTypeBoundarySensors);
	
	if(result != NULL)
	{
		result->landmark = true;
		result->broadcastCovered = 1;
		++Number_of_Landmark;
		// printf("The %dth landmark: sensor %d\n", Number_of_Landmark, result->index);
	}

	return result;
}

Sensor* Field::ReselectingFloodingSensor()
{
	unsigned int i;
	vector<Sensor*> RandomSelectingSet;
	RandomSelectingSet.clear();
	for (i = 0; i < BoundarySensorSet.size(); ++i)
	{
		Sensor *s = BoundarySensorSet[i];
		if(!s->inactive && !s->landmark && s->Rou_p == 1)
		{
			RandomSelectingSet.push_back(s);
		}
	}
	UniformRandVar urand(0, RandomSelectingSet.size()-0.1);
	double maxRou_p = 0;
	Sensor *result = NULL;
	for(i = 0; i < BoundarySensorSet.size(); ++i)
	{
		if(!BoundarySensorSet[i]->inactive && !BoundarySensorSet[i]->landmark)
		{
			if(BoundarySensorSet[i]->Rou_p > maxRou_p)
			{
				maxRou_p = BoundarySensorSet[i]->Rou_p;
				result = BoundarySensorSet[i];
			}
			if(result != NULL && result->concave)
			{
				break;
			}
		}
	}
	if(maxRou_p == 1)
	{
		result = RandomSelectingSet[urand.value()];
	}
	if(result != NULL)
	{
		result->landmark = true;
		result->broadcastCovered = 1;
		++Number_of_Landmark;
		printf("The %dth landmark: sensor %d\n", Number_of_Landmark, result->index);
	}

	return result;
}

void Field::GenerateVoronoiCell()
{
	vector<Sensor*> iTypeBoundarySensors;
	//int numOfVCELL;
	// set<Sensor*> vCellBoundary;
	const int VCELL_OFFSET = 0;
	// vector<Sensor*> 
	for (int iType = 1; iType <= MAX_BOUNDARY_TYPE; ++ iType)
	{
		cout<<"Level "<<iType<<"..."<<endl;
		// Get the ith type boundary nodes.
		for (int iSensor = 0; iSensor < nSensors; ++ iSensor)
		{
			Sensor* s = &sensorPool[iSensor];
			if ( s->isBelongToType(iType) )
			{
				iTypeBoundarySensors.push_back(s);
			}
		}
		// Generate voronoi cell.
		for (auto it = iTypeBoundarySensors.begin(); it != iTypeBoundarySensors.end(); ++ it)
		{
			Sensor* sBound = *it;
			if (!sBound->iType2Landmark[iType])
			{
				sBound->constructVoroniTile(iType, VCELL_OFFSET);
				// if come across at least 2 landmarks 
				if (sBound->iType2VoronoiCellID[iType].empty())
				{
					//++ numOfVCELL;
					// sbound is a boundary sensor of the v cells.
					cerr<<"A non-landmark node should belong to a cell."<<endl;
					cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<endl;
					exit(EXIT_FAILURE);
				}
				else
				{
					auto itFirst = sBound->iType2VoronoiCellID[iType].begin();
					sBound->iType2LandmarkIndex[iType] = *itFirst;
				}
				
			}
		}
		
		// Get the relation ship among the landmarks.
		for (auto it = iTypeBoundarySensors.begin(); it != iTypeBoundarySensors.end(); ++ it)
		{
			Sensor* sBound = *it;
			Sensor* LM_sBound    = &sensorPool[sBound->iType2LandmarkIndex[iType]];

			if (! sBound->iType2Landmark[iType] )
			{
				for (auto nbIt = sBound->neighbors.begin(); nbIt != sBound->neighbors.end(); ++ nbIt)
				{
					Sensor* nb_sBound = *nbIt;
					Sensor* LM_nb_sBound = &sensorPool[nb_sBound->iType2LandmarkIndex[iType]];
					
					// Both nb and I are on the same boundary with different landmark index
					if (!nb_sBound->iType2Landmark[iType]
					&&	nb_sBound->isBelongToType(iType) && sBound->isBelongToType(iType)
					&&  LM_nb_sBound->index           != LM_sBound->index)
					{
						LM_sBound->iTypeLandmarkNeighbors[iType].insert(LM_nb_sBound);
						LM_nb_sBound->iTypeLandmarkNeighbors[iType].insert(LM_sBound);
					}
				}
			}
		}
		//int i = 0;
		//set<int>::iterator itVID;
		//for(auto it = vCellBoundary.begin(); it != vCellBoundary.end(); ++it)
		//{
		//	Sensor *s = *it;
		//	int SIZE = s->VoronoiCellID.size();

		//	// cerr<<SIZE<<" |";

		//	int ID[128] = {0};

		//	for(i = 0, itVID = s->VoronoiCellID.begin(); itVID != s->VoronoiCellID.end(); ++itVID, ++i)
		//	{
		//		ID[i] = *itVID;
		//	}
		//	for(int i = 0; i < SIZE; ++i)
		//	{
		//		for(int j = 0; j < SIZE; ++j)
		//		{
		//			if(ID[i] != ID[j])
		//			{
		//				sensorPool[ID[i]].LandmarkNeighbors.insert(&sensorPool[ID[j]]);
		//			}
		//		}
		//	}
		//}

		// vCellBoundary.clear();
		iTypeBoundarySensors.clear();
	}

	/************************************************************************/
	/* Testing landmarks to see if they are connected.                      */
	/************************************************************************/
	cerr<<"Testing if landmarks are connected..."<<endl;

	for (int iType = 1; iType <= MAX_BOUNDARY_TYPE; ++ iType)
	{
		bool isLandmarkConnected = true;
		for (int i = 0; i < nSensors; ++ i)
		{
			Sensor* s = &sensorPool[i];
			if (s->iType2Landmark[iType] && s->iTypeLandmarkNeighbors[iType].empty())
			{
				isLandmarkConnected = false;
				break;
			}
		}
		if (!isLandmarkConnected)
		{
			cerr<<iType<<" "<<"level Landmarks are not connected!"<<endl;
			cerr<< "Error in " << __FILE__ << " at " << __LINE__ << ".\n";
			exit(EXIT_FAILURE);
		}
		else
			cerr<<iType<<" level Landmarks connected.:-D"<<endl;

		/************************************************************************/
		/*            A Good landmark should have more than 1 neighbors         */
		/************************************************************************/
		for (int i = 0; i < nSensors; ++ i)
		{
			Sensor* s = &sensorPool[i];
			if (s->iType2Landmark[iType] && s->iTypeLandmarkNeighbors[iType].size() == 1)
			{
				cout<<"****************************************************"<<endl;
				cout<<s->index<<" only have 1 neighbor..."<<endl;
				cout<<"****************************************************"<<endl;
			}
		}

		itypeAttachEdgesToEachLandmark(iType);
		itypeAttachTriangleToEachLandmark(iType);
		itypeAttachEdgesToTriangles(iType);
	}
	
	// attachEdgesToEachLandmark();
	//attachTriangleToEachLandmark();
	//attachEdgesToTriangles();

}

int Field::findEdgeInItypeNoCrossPool(const int itype, Sensor* start_, Sensor* end_)
{
	//assert(start_!=NULL && end_!=NULL);
	if (start_ == NULL || end_ == NULL)
	{
		cerr<<"The ends vertex of the lineSegment should not be null."<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	if (start_ == end_  )
	{
		cerr<<"The ends vertex of the lineSegment should not be the same."<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	LineSegment L_(-1,start_,end_);
	int idxE = -1;
	for (int iE = 0; iE < iTypeNoCrossEdgePool[itype].size(); ++iE )
	{
		if (L_ == iTypeNoCrossEdgePool[itype][iE])
		{
			idxE = iE;
			break;
		}
	}

	return idxE;

}
int Field::findEdgeInNoCrossPool(Sensor* start_, Sensor* end_)
{
	//assert(start_!=NULL && end_!=NULL);
	if (start_ == NULL || end_ == NULL)
	{
		cerr<<"The ends vertex of the lineSegment should not be null."<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	if (start_ == end_  )
	{
		cerr<<"The ends vertex of the lineSegment should not be the same."<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	LineSegment L_(-1,start_,end_);
	int idxE = -1;
	for (int iE = 0; iE < noCrossEdgePool.size(); ++iE )
	{
		if (L_ == noCrossEdgePool[iE])
		{
			idxE = iE;
			break;
		}
	}

	return idxE;

}
void Field::updateTriInfo()
{
	// the neighbors of triangle, the field of the triangle.
	for (int iType = 1; iType <= MAX_BOUNDARY_TYPE; ++ iType)
	{
		for (int i = 0; i < iTypeTriangleSet[iType].size(); ++ i)
		{
			iTypeTriangleSet[iType][i].Tri_field = this;
			iTypeTriangleSet[iType][i].iTypeEstablishTriNeighborhood(iType);
		}
	}
}
LineSegment* Field::findEdgePointer(Sensor* start_, Sensor* end_)
{
	//assert(start_!=NULL && end_!=NULL);
	if (start_ == NULL || end_ == NULL )
	{
		cerr<<"The ends vertex of the lineSegment should not be null."<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	if (start_ == end_  )
	{
		cerr<<"The ends vertex of the lineSegment should not be the same."<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	LineSegment L_(-1,start_,end_);
	int idxE = -1;
	for (int iE = 0; iE < EdgeSet.size(); ++iE )
	{
		if (L_ == EdgeSet[iE])
		{
			idxE = iE;
			break;
		}
	}

	return &(EdgeSet[idxE]);


}
LineSegment* Field::findItypeEdgePointer(const int itype,Sensor* start_, Sensor* end_)
{
	//assert(start_!=NULL && end_!=NULL);
	if (start_ == NULL || end_ == NULL )
	{
		cerr<<"The ends vertex of the lineSegment should not be null."<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	if (start_ == end_  )
	{
		cerr<<"The ends vertex of the lineSegment should not be the same."<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	LineSegment L_(-1,start_,end_);
	int idxE = -1;
	for (int iE = 0; iE < iTypeEdgesSet[itype].size(); ++iE )
	{
		if (L_ ==  iTypeEdgesSet[itype][iE])
		{
			idxE = iE;
			break;
		}
	}

	return &( iTypeEdgesSet[itype][idxE]);


}
int Field::findItypeEdge(const int itype,Sensor* start_, Sensor* end_)
{
	// Return the index in the itypeEdgeset
	//assert(start_!=NULL && end_!=NULL);
	if (start_ == NULL || end_ == NULL )
	{
		cerr<<"The ends vertex of the lineSegment should not be null."<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	if (start_ == end_  )
	{
		cerr<<"The ends vertex of the lineSegment should not be the same."<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	LineSegment L_(-1,start_,end_);
	int idxE = -1;
	for (int iE = 0; iE < iTypeEdgesSet[itype].size(); ++iE )
	{
		if (L_ == iTypeEdgesSet[itype][iE])
		{
			idxE = iE;
			break;
		}
	}

	return idxE;

}
int Field::findEdge(Sensor* start_, Sensor* end_)
{

	//assert(start_!=NULL && end_!=NULL);
	if (start_ == NULL || end_ == NULL )
	{
		cerr<<"The ends vertex of the lineSegment should not be null."<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	if (start_ == end_  )
	{
		cerr<<"The ends vertex of the lineSegment should not be the same."<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	LineSegment L_(-1,start_,end_);
	int idxE = -1;
	for (int iE = 0; iE < EdgeSet.size(); ++iE )
	{
		if (L_ == EdgeSet[iE])
		{
			idxE = iE;
			break;
		}
	}

	return idxE;

}
Triangle* Field::findItypeTrianglePointer(const int itype,Sensor* t1_, Sensor* t2_, Sensor* t3_)
{

	if (t1_ == NULL || t2_ == NULL || t3_ == NULL)
	{
		cerr<<"The vertex of the triangle should not be null."<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	Triangle T_(-1,t1_,t2_,t3_);
	int idxT = -1;
	for (int iT = 0; iT < iTypeTriangleSet[itype].size(); ++ iT)
	{
		if (T_ == iTypeTriangleSet[itype][iT])
		{
			idxT = iT;
			break;
		}
	}

	return &(iTypeTriangleSet[itype][idxT]);


}
Triangle* Field::findTrianglePointer(Sensor* t1_, Sensor* t2_, Sensor* t3_)
{
	if (t1_ == NULL || t2_ == NULL || t3_ == NULL)
	{
		cerr<<"The vertex of the triangle should not be null."<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	Triangle T_(-1,t1_,t2_,t3_);
	int idxT = -1;
	for (int iT = 0; iT < TriangleSet.size(); ++ iT)
	{
		if (T_ == TriangleSet[iT])
		{
			idxT = iT;
			break;
		}
	}

	return &(TriangleSet[idxT]);

}

int Field::findItypeTriangle(const int itype,Sensor* t1_, Sensor* t2_, Sensor* t3_)
{
	if (t1_ == NULL || t2_ == NULL || t3_ == NULL)
	{
		cerr<<"The vertex of the triangle should not be null."<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	Triangle T_(-1,t1_,t2_,t3_);
	int idxT = -1;
	for (int iT = 0; iT < iTypeTriangleSet[itype].size(); ++ iT)
	{
		if (T_ == iTypeTriangleSet[itype][iT])
		{
			idxT = iT;
			break;
		}
	}

	return idxT;
}
int Field::findTriangle(Sensor* t1_, Sensor* t2_, Sensor* t3_)
{
	if (t1_ == NULL || t2_ == NULL || t3_ == NULL)
	{
		cerr<<"The vertex of the triangle should not be null."<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	Triangle T_(-1,t1_,t2_,t3_);
	int idxT = -1;
	for (int iT = 0; iT < TriangleSet.size(); ++ iT)
	{
		if (T_ == TriangleSet[iT])
		{
			idxT = iT;
			break;
		}
	}

	return idxT;

}
void Field::itypeAttachEdgesToEachLandmark(const int itype)
{
	for (int i=0; i<nSensors;++i)
	{
		Sensor* s = &(sensorPool[i]);
		s->iTypeAttachedEdges[itype].clear();
	}
	/* attach edges to the node                                             */
	/************************************************************************/
	for (int i_edge = 0; i_edge < iTypeEdgesSet[itype].size();++i_edge)
	{
		Sensor* startV = iTypeEdgesSet[itype][i_edge].start;
		Sensor* endV = iTypeEdgesSet[itype][i_edge].end;

		assert(startV->iType2Landmark[itype] && endV->iType2Landmark[itype]);

		startV->iTypeAttachedEdges[itype].insert(&iTypeEdgesSet[itype][i_edge]);
		endV  ->iTypeAttachedEdges[itype].insert(&iTypeEdgesSet[itype][i_edge]);

	}

}
void Field::attachEdgesToEachLandmark()
{

	for (int i=0; i<nSensors;++i)
	{
		Sensor* s = &(sensorPool[i]);
		s->attachedEdges.clear();
	}
	/* attach edges to the node                                             */
	/************************************************************************/
	for (int i_edge = 0; i_edge < EdgeSet.size();++i_edge)
	{
		Sensor* startV = EdgeSet[i_edge].start;
		Sensor* endV = EdgeSet[i_edge].end;
		assert(startV->landmark && endV->landmark);

		startV->attachedEdges.insert(&EdgeSet[i_edge]);
		endV->attachedEdges.insert(&EdgeSet[i_edge]);

	}

}
void Field::itypeAttachattachEachSensorWithSickLandmarks(const int itype)
{
	for (int i = 0; i< iType2LandmarkSet[itype].size(); ++i)
	{
		Sensor* lm = iType2LandmarkSet[itype][i];
		lm->iTypeSickLandmarks[itype].clear();
	}

	itypeAttachEdgesToEachLandmark(itype);

	for (int i = 0; i< iType2LandmarkSet[itype].size(); ++i)
	{
		Sensor* lm = iType2LandmarkSet[itype][i];

		for (auto itEd = lm->iTypeAttachedEdges[itype].begin(); itEd!=lm->iTypeAttachedEdges[itype].end();
			++ itEd)
		{
			LineSegment* ed = *itEd;
			if (ed->usedtimes < 2)
			{
				lm->iTypeSickLandmarks[itype].insert(ed->start);
				lm->iTypeSickLandmarks[itype].insert(ed->end);
			}

		}

		lm->iTypeSickLandmarks[itype].erase(lm);

	}

}
void Field::attachEachSensorWithSickLandmarks()
{
	for (int i = 0; i< LandmarkSet.size(); ++i)
	{
		Sensor* lm = LandmarkSet[i];
		lm->sickLandmarks.clear();
	}

	attachEdgesToEachLandmark();

	for (int i = 0; i< LandmarkSet.size(); ++i)
	{
		Sensor* lm = LandmarkSet[i];

		for (auto itEd = lm->attachedEdges.begin(); itEd!=lm->attachedEdges.end();
			++ itEd)
		{
			LineSegment* ed = *itEd;
			if (ed->usedtimes < 2)
			{
				lm->sickLandmarks.insert(ed->start);
				lm->sickLandmarks.insert(ed->end);
			}
			
		}

		lm->sickLandmarks.erase(lm);

	}


	
}

void Field::buildEdgeAttachedLandmarkNb()
{
	for (int i = 0; i< LandmarkSet.size(); ++i)
	{
		Sensor* lm = LandmarkSet[i];
		lm->edgeAttachedLandmarkNb.clear();
	}

	attachEdgesToEachLandmark();

	for (int i = 0; i< LandmarkSet.size(); ++i)
	{
		Sensor* lm = LandmarkSet[i];

		for (auto itEd = lm->attachedEdges.begin(); itEd!=lm->attachedEdges.end();
			++ itEd)
		{
			LineSegment* ed = *itEd;
			
			lm->edgeAttachedLandmarkNb.insert(ed->start);
			lm->edgeAttachedLandmarkNb.insert(ed->end);

		}

		lm->edgeAttachedLandmarkNb.erase(lm);

	}

}
void Field::itypeAttachTriangleToEachLandmark(const int itype)
{
	for (int i=0; i<nSensors; ++i)
	{
		Sensor* s = &(sensorPool[i]);
		s->iTypeAttachedTriangles[itype].clear();
	}
	for (int iT = 0; iT < iTypeTriangleSet[itype].size(); ++iT)
	{
		Sensor* v1 = iTypeTriangleSet[itype][iT].t1;
		Sensor* v2 = iTypeTriangleSet[itype][iT].t2;
		Sensor* v3 = iTypeTriangleSet[itype][iT].t3;

		assert(v1->iType2Landmark[itype] && v2->iType2Landmark[itype] && v3->iType2Landmark[itype]);

		v1->iTypeAttachedTriangles[itype].insert(&iTypeTriangleSet[itype][iT]);
		v2->iTypeAttachedTriangles[itype].insert(&iTypeTriangleSet[itype][iT]);
		v3->iTypeAttachedTriangles[itype].insert(&iTypeTriangleSet[itype][iT]);

	}

}
void Field::attachTriangleToEachLandmark()
{
	for (int i=0; i<nSensors; ++i)
	{
		Sensor* s = &(sensorPool[i]);
		s->attachedTriangles.clear();
	}
	for (int iT = 0; iT < TriangleSet.size(); ++iT)
	{
		Sensor* v1 = TriangleSet[iT].t1;
		Sensor* v2 = TriangleSet[iT].t2;
		Sensor* v3 = TriangleSet[iT].t3;

		assert(v1->landmark && v2->landmark && v3->landmark);

		v1->attachedTriangles.insert(&TriangleSet[iT]);
		v2->attachedTriangles.insert(&TriangleSet[iT]);
		v3->attachedTriangles.insert(&TriangleSet[iT]);

	}
}
void Field::itypeAttachEdgesToTriangles(const int itype)
{// to change....
	// a edge is attached with triangles.
	// EdgeSet and TriangleSet are used.
	for (int iE = 0; iE < iTypeEdgesSet[itype].size();++iE)
	{
		iTypeEdgesSet[itype][iE].itypeEdgeAttachedTriangles[itype].clear();
	}
	for (int iT = 0; iT < iTypeTriangleSet[itype].size();++iT)
	{
		Sensor* v1 = iTypeTriangleSet[itype][iT].t1;
		Sensor* v2 = iTypeTriangleSet[itype][iT].t2;
		Sensor* v3 = iTypeTriangleSet[itype][iT].t3;

		int e12 = findItypeEdge(itype,v1,v2);
		int e23 = findItypeEdge(itype,v2,v3);
		int e13 = findItypeEdge(itype,v1,v3);

		assert(e12>-1&&e23>-1&&e13>-1);//edges must exist.

		iTypeEdgesSet[itype][e12].itypeEdgeAttachedTriangles[itype].insert(&iTypeTriangleSet[itype][iT]);
		iTypeEdgesSet[itype][e23].itypeEdgeAttachedTriangles[itype].insert(&iTypeTriangleSet[itype][iT]);
		iTypeEdgesSet[itype][e13].itypeEdgeAttachedTriangles[itype].insert(&iTypeTriangleSet[itype][iT]);
	}
}

void Field::attachEdgesToTriangles()
{
	// a edge is attached with triangles.
	// EdgeSet and TriangleSet are used.
	for (int iE = 0; iE < EdgeSet.size();++iE)
	{
		EdgeSet[iE].attachedTriangles.clear();
	}
	for (int iT = 0; iT < TriangleSet.size();++iT)
	{
		Sensor* v1 = TriangleSet[iT].t1;
		Sensor* v2 = TriangleSet[iT].t2;
		Sensor* v3 = TriangleSet[iT].t3;

		int e12 = findEdge(v1,v2);
		int e23 = findEdge(v2,v3);
		int e13 = findEdge(v1,v3);

		assert(e12>-1&&e23>-1&&e13>-1);//edges must exist.

		EdgeSet[e12].attachedTriangles.insert(&TriangleSet[iT]);
		EdgeSet[e23].attachedTriangles.insert(&TriangleSet[iT]);
		EdgeSet[e13].attachedTriangles.insert(&TriangleSet[iT]);
	}
}
bool Field::isLeadingToIntersectionItype(const int itype, LineSegment& tmpLine)
{
	// the path on the boundary should not be empty.
	const int minPathLength = 2;
	assert(tmpLine.itypeBoundaryPath[itype].size() > minPathLength);
	assert(tmpLine.start->iType2Landmark[itype] && tmpLine.end->iType2Landmark[itype]);

	bool intersect;
	int posInPool = findEdgeInItypeNoCrossPool(itype, tmpLine.start,tmpLine.end);
	if (posInPool == -1)
	{
		intersect = true;
	}
	else
		intersect = false;

	return intersect;

}
bool Field::isLeadingToIntersection(LineSegment& tmpLine)
{

	// the path on the boundary should not be empty.
	const int minPathLength = 2;
	assert(tmpLine.boundaryPath.size() > minPathLength);
	assert(tmpLine.start->landmark && tmpLine.end->landmark);

	bool intersect;
	int posInPool = findEdgeInNoCrossPool(tmpLine.start,tmpLine.end);
	if (posInPool == -1)
	{
		intersect = true;
	}
	else
		intersect = false;

	return intersect;

}
void Field::egdesOfLandmarks()
{
	cerr<<"Building no crossing edge pool..."<<endl;
	buildNoCrossingEdgePool();
	cerr<<"No crossing edge pool is built."<<endl;
	cerr<<"Partial virtual edges extraction..."<<endl;

	for (int iType = 1; iType <= MAX_BOUNDARY_TYPE; ++ iType)
	{
		// Connect the landmarks.
		int lineIdx = 0;

		for (int iLm = 0; iLm < iType2LandmarkSet[iType].size(); ++ iLm)
		{
			Sensor* lmS = iType2LandmarkSet[iType][iLm];
			set<Sensor*>noIntersectionLM_Nb;
			if (lmS->iTypeAttachedEdges[iType].empty())
			{
				for (auto it_lmS_nb = lmS->iTypeLandmarkNeighbors[iType].begin(); 
					it_lmS_nb != lmS->iTypeLandmarkNeighbors[iType].end();++ it_lmS_nb)
				{
					LineSegment LTmp(-1,lmS,*it_lmS_nb);
					LTmp.floodStart2EndPathOnItypeBoundary(iType);
					if (!isLeadingToIntersectionItype(iType,LTmp) && (*it_lmS_nb)->iTypeAttachedEdges[iType].empty())
					{
						noIntersectionLM_Nb.insert(*it_lmS_nb);
					}

				}
				// find the landmark with greatest importance.
				if (!noIntersectionLM_Nb.empty())
				{
					auto it_GreatestImportanceLM_Nb = max_element(noIntersectionLM_Nb.begin(),noIntersectionLM_Nb.end(),[&](Sensor* s1, Sensor* s2) -> bool
					{
						return (s1->Rou_p < s2->Rou_p);
					});

					LineSegment goodLine(lineIdx,lmS,*it_GreatestImportanceLM_Nb);
					goodLine.floodStart2EndPathOnItypeBoundary(iType);
					iTypeEdgesSet[iType].push_back(goodLine);
					++ lineIdx;
					itypeAttachEdgesToEachLandmark(iType);
				}
			}
		}
		cerr<<iTypeEdgesSet[iType].size()<<" virtual edges extracted in level "<<iType<<"."<<endl;
	}
	
}

bool Field::isEveryItypeEdgeSaturated(const int itype)
{
	bool saturated = true;
	if (iTypeEdgesSet[itype].empty())
	{
		cerr<<"No partial virtual edges constructed"<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	else
	{
		for (int iE = 0; iE < iTypeEdgesSet[itype].size(); ++iE)
		{
			if (iTypeEdgesSet[itype][iE].usedtimes < 2)
			{
				saturated = false;
				break;
			}
		}

	}

	return saturated;
}
bool Field::isEveryEdgeSaturated()
{
	bool saturated = true;
	if (EdgeSet.empty())
	{
		cerr<<"No partial virtual edges constructed"<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	else
	{
		for (int iE = 0; iE < EdgeSet.size(); ++iE)
		{
			if (EdgeSet[iE].usedtimes < 2)
			{
				saturated = false;
				break;
			}
		}

	}

	return saturated;

}

void Field::testingTriangleOutPut()
{
	for (unsigned int i = 0; i < TriangleSet.size(); ++i)
	{
		TriangleSet[i].GetNormalVector();
		TriangleSet[i].Tri_field = this;
		TriangleSet[i].EstablishNeighborhood();
	}

// 	int wrongEdgeNum = 0;
// 	int needtoprocess = 0;
// 	for(unsigned int i  = 0; i < EdgeSet.size(); ++i)
// 	{
// 		if(EdgeSet[i].usedtimes > 2)
// 			++wrongEdgeNum;
// 		if(EdgeSet[i].usedtimes < 2)
// 			++needtoprocess;
// 	}
// 	printf("The size of edge set is %d\n", EdgeSet.size());
// 	printf("The size of refining edge set is %d\n", RefineEdgeSet.size());
// 	printf("The number of wrong edge is %d\n", wrongEdgeNum);

	
	//Computing distance error
	//set<Sensor*> OuterBoundaryNodes;
	double sum_dist_error = 0;
	//refine triangles' normal vector
	for(unsigned int i = 0; i < TriangleSet.size(); ++i)
	{
		Triangle *t = &TriangleSet[i];
		Point reference_point;
		reference_point.x = (t->t1->location.x+t->t2->location.x+t->t3->location.x)/3;
		reference_point.y = (t->t1->location.y+t->t2->location.y+t->t3->location.y)/3;
		reference_point.z = (t->t1->location.z+t->t2->location.z+t->t3->location.z)/3;

		//confirm the notation of normal
		int xrange = 200;
		int yrange = 200;
		int zrange = 200;
		Point Orientation_center;
		int Length = LENGTH;// This is a ????
		Orientation_center.x = Length*t->Normalx+reference_point.x;
		Orientation_center.y = Length*t->Normaly+reference_point.y;
		Orientation_center.z = Length*t->Normalz+reference_point.z;
		/*
		//for S topology
		if((Orientation_center.y>=0 && Orientation_center.y <= yrange/5 && Orientation_center.z>=0 && Orientation_center.z<=zrange && Orientation_center.x>=0 && Orientation_center.x<=xrange) 
			|| (Orientation_center.y>=yrange/5*2 && Orientation_center.y<=yrange/5*3 && Orientation_center.z>=0 && Orientation_center.z<=zrange && Orientation_center.x>=0 && Orientation_center.x<=xrange) 
			|| (Orientation_center.y>=yrange/5*4 && Orientation_center.y<=yrange && Orientation_center.z>=0 && Orientation_center.z<=zrange && Orientation_center.x>=0 && Orientation_center.x<=xrange)
			|| (Orientation_center.y>=yrange/5 && Orientation_center.y<=yrange/5*2 && Orientation_center.x>=xrange/5*3 && Orientation_center.x<=xrange && Orientation_center.z>=0 && Orientation_center.z<=zrange) 
			|| (Orientation_center.x>=0 && Orientation_center.x<=xrange/5*2 && Orientation_center.y>=yrange/5*3 && Orientation_center.y<=yrange/5*4 && Orientation_center.z>=0 && Orientation_center.z<=zrange))
		{
			t->Normalx = 0-t->Normalx;
			t->Normaly = 0-t->Normaly;
			t->Normalz = 0-t->Normalz;
		}
		*/
		//for chicago airport & cube with H

		if(WithinTopology(Orientation_center))
		{
			t->Normalx = 0-t->Normalx;
			t->Normaly = 0-t->Normaly;
			t->Normalz = 0-t->Normalz;
		}
	}
	//refine again
	for(unsigned int i = 0; i < TriangleSet.size(); ++i)
	{
		Triangle *t = &TriangleSet[i];
		int matchtimes = 0;
		if(t->TriangleNeighbors.size() > 3)
		{
			matchtimes -= (t->TriangleNeighbors.size()-3);
		}
		for(set<Triangle*>::iterator TNpos = t->TriangleNeighbors.begin(); TNpos != t->TriangleNeighbors.end(); ++TNpos)
		{
			Triangle *tn = *TNpos;
			double product = t->Normalx*tn->Normalx + t->Normaly*tn->Normaly + t->Normalz*tn->Normalz;
			if(product < 0)
			{
				++matchtimes;
			}
		}
		if(matchtimes > 2)
		{
			t->Normalx = 0-t->Normalx;
			t->Normaly = 0-t->Normaly;
			t->Normalz = 0-t->Normalz;
		}
	}
	//////////////////////////////////////////////////////////////////////////
	for(int i = 0; i < nSensors; ++i)
	{
		if(!sensorPool[i].landmark /*&& sensorPool[i].onEdge_Degree != 0*/)
		{
			Sensor *s = &sensorPool[i];
			double sum_distance = 0xffff;
			Triangle *t = NULL;
			for(unsigned int j = 0; j < TriangleSet.size(); ++j)
			{
				double sd = 0;
				sd += s->location.distance(TriangleSet[j].t1->location);
				sd += s->location.distance(TriangleSet[j].t2->location);
				sd += s->location.distance(TriangleSet[j].t3->location);
				if(sd < sum_distance)
				{
					sum_distance = sd;
					t = &TriangleSet[j];
				}
			}
			

			Point reference_point;
			reference_point.x = (t->t1->location.x+t->t2->location.x+t->t3->location.x)/3;
			reference_point.y = (t->t1->location.y+t->t2->location.y+t->t3->location.y)/3;
			reference_point.z = (t->t1->location.z+t->t2->location.z+t->t3->location.z)/3;

			double ve_x = s->location.x-reference_point.x;
			double ve_y = s->location.y-reference_point.y;
			double ve_z = s->location.z-reference_point.z;
			double ve_product = ve_x*t->Normalx + ve_y*t->Normaly + ve_z*t->Normalz;
			if(ve_product > 0)
			{
				s->Out = true;
				++num_of_outer_boundary_surface;
				double A = fabs(t->Normalx*(s->location.x-t->t1->location.x) + t->Normaly*(s->location.y - t->t1->location.y) + t->Normalz*(s->location.z - t->t1->location.z));
				double B = sqrt(t->Normalx*t->Normalx+t->Normaly*t->Normaly+t->Normalz*t->Normalz);
				double d = A/B;
				sum_dist_error += d;
			}
		}
	}
	Avg_distance_err = sum_dist_error/nSensors;
	cerr<<"Triangle set size is "<<TriangleSet.size()<<endl;
	dumpTriangle("triangle");
	dumpTriangleWithID("triangle.index");
	dumpEdges("EdgeSetedges");

}

void Field::buildNoCrossingEdgePool()
{
	// Find all no crossing edges.
	// Connect the landmarks.
	int lineIdx = 0;
	vector<Sensor*> iTypeSensors;
	vector<Sensor*> iTypeLandmarks;
	// vector<LineSegment> LandmarkEdges;
	for (int iType = 1; iType <= MAX_BOUNDARY_TYPE; ++ iType)
	{
		// Get the itype boundary nodes.
		for (int iSensor = 0; iSensor < nSensors; ++ iSensor)
		{
			Sensor* s = &sensorPool[iSensor];
			if (s->isBelongToType(iType))
			{
				iTypeSensors.push_back(s);
			}
		}
		// Find the landmarks in itype sensors.
		for (auto it = iTypeSensors.begin(); it != iTypeSensors.end(); ++ it)
		{
			Sensor* sbound = *it;
			if (sbound->iType2Landmark[iType])
			{
				iTypeLandmarks.push_back(sbound);
			}
		}
		// Connecting the landmarks.
		for (auto itLM = iTypeLandmarks.begin(); itLM != iTypeLandmarks.end(); ++ itLM)
		{
			Sensor* sLM = *itLM;
			set<Sensor*>::iterator it_sLM_nb = sLM->iTypeLandmarkNeighbors[iType].begin();
			// bool isGoodPath;

			while (it_sLM_nb != sLM->iTypeLandmarkNeighbors[iType].end())
			{
				Sensor* sLM_nb = *it_sLM_nb;
				// nb should be a landmark.
				assert(sLM_nb->iType2Landmark[iType]);

				/*******************************************/
				// very important, for increment and erase.				
				set<Sensor*>::iterator current_it_sLM_nb = it_sLM_nb ++;
				/*******************************************/
			
				// a good path is a good edge and should be saved.

				LineSegment LM_line(lineIdx,sLM,sLM_nb);
				auto it = find_if(iTypeNoCrossEdgePool[iType].begin(),iTypeNoCrossEdgePool[iType].end(),[&](LineSegment L)
				{
					return (L == LM_line);
				});
				// can not find the edge LM_line
				if (it == iTypeNoCrossEdgePool[iType].end())
				{
					LM_line.line_Field = this;
					LM_line.floodStart2EndPathOnItypeBoundary(iType);
					cerr<<LM_line.itypeBoundaryPath[iType].size()<<" | ";
					LM_line.valid = true;
					iTypeNoCrossEdgePool[iType].push_back(LM_line);
					
					++ lineIdx;
				}	
			}	
		}

		iTypeLandmarks.clear();
		iTypeSensors.clear();
		
		cerr<<endl;
		cerr<<iTypeNoCrossEdgePool[iType].size()<<" Edges are found in level "<<iType<<"."<<endl;
		

	}

	cout<<"Removing crossing edges..."<<endl;

	// Till now, all the legal edges are found, and the edges should not cross...
	for (int iType = 1; iType <= MAX_BOUNDARY_TYPE; ++ iType)
	{
		for (int iE = 0; iE < iTypeNoCrossEdgePool[iType].size(); ++ iE)
		{
			for ( int jE = iE + 1; jE < iTypeNoCrossEdgePool[iType].size(); ++ jE )
			{
				if (iTypeNoCrossEdgePool[iType][iE].TheSameLine(&(iTypeNoCrossEdgePool[iType][jE])))
				{
					// the two line should not be the same line
					cerr<<"Same lines in the EdgeSet."<<endl;
					cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<endl;
					exit(EXIT_FAILURE);
				}
				else
				{
					// the two lines should not sharing a vertex.
					/*LineSegment IE_L(EdgeSet[iE]);
					LineSegment JE_L(EdgeSet[jE]);*/ //the two lines are for debug

					if (iTypeNoCrossEdgePool[iType][iE].valid && iTypeNoCrossEdgePool[iType][jE].valid 
						&& !iTypeNoCrossEdgePool[iType][iE].isSharingAVertexWith(iTypeNoCrossEdgePool[iType][jE]))
					{
						bool crossing = false;

						// Be careful here, i and j is different.
						// I think the crossing 

						assert(iTypeNoCrossEdgePool[iType][iE].itypeBoundaryPath[iType].size() > 1);

						for (auto it_s_iEPath = iTypeNoCrossEdgePool[iType][iE].itypeBoundaryPath[iType].begin();
							it_s_iEPath != iTypeNoCrossEdgePool[iType][iE].itypeBoundaryPath[iType].end()-1;++ it_s_iEPath)
						{
							Sensor* s_iEPath = *it_s_iEPath;

							assert(iTypeNoCrossEdgePool[iType][jE].itypeBoundaryPath[iType].size() > 1);

							for (auto it_s_jEPath = iTypeNoCrossEdgePool[iType][jE].itypeBoundaryPath[iType].begin();
								it_s_jEPath != iTypeNoCrossEdgePool[iType][jE].itypeBoundaryPath[iType].end()-1;++ it_s_jEPath)
							{

								Sensor* s_jEpath = *it_s_jEPath;

								// the two sensors should not be boundary neighbors
								auto it_nb = find_if(s_iEPath->neighbors.begin(),s_iEPath->neighbors.end(),[&](Sensor* s) -> bool
								{
									return(s == s_jEpath);
								});

								if (s_jEpath == s_iEPath || it_nb != s_iEPath->neighbors.end())
								{

									crossing = true;// For breaking...
									// crossing process
									//1. set the iE edge invalid
									//2. cancel the landmark neighbors of iE start and end.
									//break the for loops that processing path..
									Sensor* sIE_Start = iTypeNoCrossEdgePool[iType][iE].start;
									Sensor* sIE_end   = iTypeNoCrossEdgePool[iType][iE].end;

									Sensor* sJE_Start = iTypeNoCrossEdgePool[iType][jE].start;
									Sensor* sJE_end   = iTypeNoCrossEdgePool[iType][jE].end;

									if (sIE_Start->iTypeLandmarkNeighbors[iType].size() >2 
										&& sIE_end->iTypeLandmarkNeighbors[iType].size() >2)
									{
										iTypeNoCrossEdgePool[iType][iE].valid = false;

										//And they should not be landmark neighbors.
										sIE_Start->iTypeLandmarkNeighbors[iType].erase(sIE_end);
										sIE_end->iTypeLandmarkNeighbors[iType].erase(sIE_Start);
									}
									else if (sJE_Start->iTypeLandmarkNeighbors[iType].size() >2 
										&&   sJE_end->iTypeLandmarkNeighbors[iType].size() >2)
									{
										iTypeNoCrossEdgePool[iType][jE].valid = false;

										//And they should not be landmark neighbors.
										sJE_Start->iTypeLandmarkNeighbors[iType].erase(sJE_end);
										sJE_end->iTypeLandmarkNeighbors[iType].erase(sJE_Start);
									}
									else
										cerr<<"Can not delete...maybe not right."<<endl;


									break;

								}

							}

							if (crossing)
							{
								//break the second path for loop.
								break;
							}

						}
					}
				}


			}
		}

		//Remove the invalid edges.
		iTypeNoCrossEdgePool[iType].erase(remove_if(iTypeNoCrossEdgePool[iType].begin(),iTypeNoCrossEdgePool[iType].end(),
			[&](LineSegment L) -> bool
		{
			return(! L.valid);
		}),iTypeNoCrossEdgePool[iType].end());

		// the edges in EdgeSet supposed to be uncrossed.
		cerr<<iTypeNoCrossEdgePool[iType].size()<<" edges left after removing the crossing edges in level "<<iType<<"."<<endl;

		/* Testing landmarks to see if they are connected.                      */
		/************************************************************************/
		cerr<<"Re Testing if landmarks are connected..."<<endl;
		bool isLandmarkConnected = true;
		for (int i = 0; i < nSensors; ++ i)
		{
			Sensor* s = &sensorPool[i];
			if (s->iType2Landmark[iType] && s->iTypeLandmarkNeighbors[iType].empty())
			{
				isLandmarkConnected = false;
				break;
			}
		}
		if (!isLandmarkConnected)
		{
			cerr<<"Landmarks are not connected!"<<endl;
			cerr<< "Error in " << __FILE__ << " at " << __LINE__ << ".\n";
			exit(EXIT_FAILURE);
		}
		else
			cerr<<"Landmarks connected.:-D"<<endl;

		/************************************************************************/
		/*            A Good landmark should have more than 1 neighbors         */
		/************************************************************************/
		for (int i = 0; i < nSensors; ++ i)
		{
			Sensor* s = &sensorPool[i];
			if (s->iType2Landmark[iType] && s->iTypeLandmarkNeighbors[iType].size() == 1)
			{
				cout<<"****************************************************"<<endl;
				cout<<s->index<<" only have 1 neighbor..."<<endl;
				cout<<"****************************************************"<<endl;
				exit(EXIT_FAILURE);
			}
		}
	}

	// A node is attached with edges.
	/************************************************************************/
	// attachEdgesToEachLandmark();
}

void Field::landmarkTriangulation()
{
	for (int iType = 1; iType <= MAX_BOUNDARY_TYPE; ++ iType)
	{
		cerr<<iType<<" level "<<"Triangles Evolution..."<<endl;
		int newLIdx = iTypeEdgesSet[iType].size()+1;
		int iE = 0;
		int triangleIdx = 0;//triangle index...
		int oldTriangleSetSize = -1;
		int newTriangleSetSize = 0;
		int roundnum = 0;
		const int limitRound = 6;
		while (!isEveryItypeEdgeSaturated(iType) && roundnum < limitRound)
		{
			roundnum ++;
			oldTriangleSetSize = iTypeTriangleSet[iType].size();
			cerr<<oldTriangleSetSize<<" triangles in the beginning of loop in level "<<iType<<" ..."<<endl;
			for (iE = 0; iE < iTypeEdgesSet[iType].size(); ++ iE)
			{
				// cerr<<iE<<" loop"<<endl;

				// 			if (iE == 110)
				// 			{
				// 				cerr<<"test"<<endl;
				// 			}

				if (iTypeEdgesSet[iType][iE].usedtimes < 2)
				{
					Sensor* startVertex = iTypeEdgesSet[iType][iE].start;
					Sensor* endVertex   = iTypeEdgesSet[iType][iE].end;



					set<Sensor*>ithEdgeNb;
					ithEdgeNb.insert(startVertex->iTypeLandmarkNeighbors[iType].begin(),
									 startVertex->iTypeLandmarkNeighbors[iType].end() );
					ithEdgeNb.insert(endVertex->iTypeLandmarkNeighbors[iType].begin(), 
									 endVertex->iTypeLandmarkNeighbors[iType].end()   );

					ithEdgeNb.erase(startVertex);
					ithEdgeNb.erase(endVertex);

					set<Sensor*>legalEdgeNb;


					for (auto itThirdS = ithEdgeNb.begin(); itThirdS != ithEdgeNb.end(); ++ itThirdS)
					{
						Sensor* thirdS = *itThirdS;

						int s2thirdPos = findItypeEdge(iType,startVertex,thirdS);
						int e2thirdPos = findItypeEdge(iType,endVertex,thirdS);

						if (s2thirdPos>-1 && iTypeEdgesSet[iType][s2thirdPos].usedtimes == 2)
						{
							continue;
						}
						if (e2thirdPos>-1 && iTypeEdgesSet[iType][e2thirdPos].usedtimes == 2)
						{
							continue;
						}

						LineSegment ls2third(-1,thirdS,startVertex);
						ls2third.floodStart2EndPathOnItypeBoundary(iType);


						if ( isLeadingToIntersectionItype(iType,ls2third))
						{
							continue;
						}

						LineSegment le2third(-1,thirdS,endVertex);
						le2third.floodStart2EndPathOnItypeBoundary(iType);

						if (isLeadingToIntersectionItype(iType,le2third))
						{
							continue;
						}

						legalEdgeNb.insert(thirdS);

					}
					if (!legalEdgeNb.empty())
					{
						auto it_greatestThird = max_element(legalEdgeNb.begin(),legalEdgeNb.end(),[&](Sensor* s1, Sensor* s2) -> bool
						{
							return (s1->Rou_p < s2->Rou_p);
						});

						Sensor* greatestThird = *it_greatestThird;
						int triPos = findItypeTriangle(iType,startVertex,endVertex,greatestThird);
						// find in the constructed set.
						Triangle LM_triangle(triangleIdx,startVertex,endVertex,greatestThird);

						auto it_removed = find_if(iTypeRemovedTris[iType].begin(),iTypeRemovedTris[iType].end(),[&](Triangle& tInS)
						{
							return (tInS == LM_triangle);
						});
						if (triPos == -1 && it_removed == iTypeRemovedTris[iType].end() )
						{

							iTypeTriangleSet[iType].push_back(LM_triangle);
							++ triangleIdx;

							// 						if ( (startVertex->index == 23913 && greatestThird->index == 28170)
							// 							||(startVertex->index == 28170 && greatestThird->index == 23913))
							// 						{
							// 							cerr<<"Problems..."<<endl;
							// 						}
							// 
							// 						if ( (endVertex->index == 23913 && greatestThird->index == 28170)
							// 							||(endVertex->index == 28170 && greatestThird->index == 23913))
							// 						{
							// 							cerr<<"Problems..."<<endl;
							// 						}

							int lstart2thirdPos = findItypeEdge(iType,startVertex,greatestThird);
							int lend2thirdPos   = findItypeEdge(iType,endVertex,greatestThird);

							if (lstart2thirdPos == -1 && lend2thirdPos > -1)
							{
								// 							int  u1 = EdgeSet[iE].usedtimes;
								// 							int  u2 = EdgeSet[lend2thirdPos].usedtimes;

								iTypeEdgesSet[iType][iE].usedtimes ++;
								iTypeEdgesSet[iType][lend2thirdPos].usedtimes++;

								LineSegment lS2Third(newLIdx,startVertex,greatestThird);
								lS2Third.floodStart2EndPathOnItypeBoundary(iType);
								lS2Third.usedtimes ++;
								iTypeEdgesSet[iType].push_back(lS2Third);
								++ newLIdx;

							}
							else if (lend2thirdPos == -1 && lstart2thirdPos > -1)
							{
								// 							int u1 = EdgeSet[iE].usedtimes;
								// 							int u2 = EdgeSet[lstart2thirdPos].usedtimes;

								iTypeEdgesSet[iType][iE].usedtimes ++;
								iTypeEdgesSet[iType][lstart2thirdPos].usedtimes++;

								LineSegment le2Third(newLIdx,endVertex,greatestThird);
								le2Third.floodStart2EndPathOnItypeBoundary(iType);
								(le2Third.usedtimes) ++;
								iTypeEdgesSet[iType].push_back(le2Third);
								++ newLIdx;

							}
							else if (lstart2thirdPos == -1 && lend2thirdPos == -1)
							{
								/*int u1 = EdgeSet[iE].usedtimes;*/

								iTypeEdgesSet[iType][iE].usedtimes ++;

								LineSegment lS2Third(newLIdx,startVertex,greatestThird);
								lS2Third.floodStart2EndPathOnItypeBoundary(iType);
								lS2Third.usedtimes ++;
								iTypeEdgesSet[iType].push_back(lS2Third);
								++ newLIdx;

								LineSegment le2Third(newLIdx,endVertex,greatestThird);
								le2Third.floodStart2EndPathOnItypeBoundary(iType);
								le2Third.usedtimes ++;
								iTypeEdgesSet[iType].push_back(le2Third);
								++ newLIdx;
							}
							else
							{
								// 							int u1 = EdgeSet[iE].usedtimes;
								// 							int u2 = EdgeSet[lstart2thirdPos].usedtimes;
								// 							int u3 = EdgeSet[lend2thirdPos].usedtimes;

								iTypeEdgesSet[iType][lstart2thirdPos].usedtimes ++;
								iTypeEdgesSet[iType][lend2thirdPos].usedtimes ++;
								iTypeEdgesSet[iType][iE].usedtimes ++;
							}

							itypeAttachEdgesToEachLandmark(iType);
							itypeAttachTriangleToEachLandmark(iType);
							itypeAttachEdgesToTriangles(iType);
							
// 							attachEdgesToEachLandmark();
// 							attachTriangleToEachLandmark();
// 							attachEdgesToTriangles();


						}

					}



				}
			}

			newTriangleSetSize = iTypeTriangleSet[iType].size();
			cerr<<newTriangleSetSize<<" triangles in the ending of loop in level "<<iType<<" ..."<<endl;
			// remove the teras...
			itypeAttachEdgesToEachLandmark(iType);
			itypeAttachTriangleToEachLandmark(iType);
			itypeAttachEdgesToTriangles(iType);

			//removeItypeTeras(iType);

			//cerr<<iTypeTriangleSet[iType].size()<<" triangles left after removing teras in level "<<iType<<" ..."<<endl;
			
// 			itypeAttachEdgesToEachLandmark(iType);
// 			itypeAttachTriangleToEachLandmark(iType);
// 			itypeAttachEdgesToTriangles(iType);
		}

		// find abnormal edges...
		int zeroTimes = 0;
		int oneTimes = 0;
		int twoTimes = 0;
		for (int iE = 0; iE < iTypeEdgesSet[iType].size(); ++ iE)
		{

			if ( iTypeEdgesSet[iType][iE].usedtimes == 0)
			{
				++ zeroTimes;
			}
			else if ( iTypeEdgesSet[iType][iE].usedtimes == 1)
			{
				++ oneTimes;
			}
			else if ( iTypeEdgesSet[iType][iE].usedtimes == 2)
			{
				++ twoTimes;
			}
			else
				cerr<<"Impossible edge used times."<<endl;
		}
		cerr<<"In level "<<iType<<"..."<<endl;
		cerr<<zeroTimes<<" edges used 0 time."<<endl;
		cerr<<oneTimes<< " edges used 1 time."<<endl;
		cerr<<twoTimes<< " edges used 2 times."<<endl;

		cerr<<test3ItypeEdges1usedtimeHole(iType)<<" 1-1-1 holes"<<endl;


		//test if the edges are in the edge pool...
		int counterAbnormalEdges = 0;
		for (int ie = 0; ie < iTypeEdgesSet[iType].size(); ++ie)
		{
			int pos = findEdgeInItypeNoCrossPool(iType,iTypeEdgesSet[iType][ie].start,iTypeEdgesSet[iType][ie].end);
			if (pos == -1)
			{
				++ counterAbnormalEdges;
			}
		}

		cout<<counterAbnormalEdges<<" edges are not in the no cross pool in level "<<iType<<" ."<<endl;

		//test 
		int counterTriAbnormalEdges = 0;
		for (int itr = 0; itr < iTypeTriangleSet[iType].size(); ++itr)
		{
			int e12 = findEdgeInItypeNoCrossPool(iType,iTypeTriangleSet[iType][itr].t1,iTypeTriangleSet[iType][itr].t2);
			int e23 = findEdgeInItypeNoCrossPool(iType,iTypeTriangleSet[iType][itr].t2,iTypeTriangleSet[iType][itr].t3);
			int e13 = findEdgeInItypeNoCrossPool(iType,iTypeTriangleSet[iType][itr].t1,iTypeTriangleSet[iType][itr].t3);

			if (e12==-1 || e23== -1 ||e13==-1)
			{
				++ counterTriAbnormalEdges;
			}
		}

		cout<<counterTriAbnormalEdges<<" triangle edges are not in the no cross pool in level "<<iType<<" ."<<endl;

		cerr<<iTypeTriangleSet[iType].size()<< " triangles in the beginning of the filling hole step in level "<<iType<<" ."<<endl;
		cerr<<iTypeEdgesSet[iType].size()<<" edges in the beginning of the filling hole step in level "<<iType<<" ."<<endl;



		// Find the holes...
		// there is no impossible holes.
		// collect the sick edges.

		int lineIndex = iTypeEdgesSet[iType].size() + 1;
		int triIndex  = iTypeTriangleSet[iType].size() + 1;

		itypeAttachattachEachSensorWithSickLandmarks(iType);// now every sensor has valid Landmark neighbors.
		itypeAttachEdgesToEachLandmark(iType);
		itypeAttachTriangleToEachLandmark(iType);
		itypeAttachEdgesToTriangles(iType);

		set<Sensor*> minLenSensors;
		while(true)
		{
			// find the landmarks that have sick landmarks.
			set<Sensor*>sickLm;
			for (int i=0; i<iType2LandmarkSet[iType].size(); ++i)
			{
				Sensor* s = iType2LandmarkSet[iType][i];
				if (! s->iTypeSickLandmarks[iType].empty())
				{
					sickLm.insert(s);
				}
			}
			// find the 1 sick landmark.
			set<Sensor*>sickLm_1nb;
			for (int i=0; i<iType2LandmarkSet[iType].size(); ++i)
			{
				Sensor* s = iType2LandmarkSet[iType][i];
				if (s->iTypeSickLandmarks[iType].size() == 1)
				{
					sickLm_1nb.insert(s);
				}
			}
			// eliminate the sickLm_1nb
			if (! sickLm_1nb.empty())
			{
				for (auto it = sickLm_1nb.begin(); it!= sickLm_1nb.end(); ++ it)
				{
					Sensor* s1 = *it;

					if (s1->iTypeSickLandmarks[iType].empty())
					{
						set<Sensor*>sickLm_copy = sickLm;
						sickLm_copy.erase(s1);
						for (auto itMin = minLenSensors.begin(); itMin != minLenSensors.end(); ++ itMin)
						{
							sickLm_copy.erase(*itMin);
						}

						auto it_min_len = min_element(sickLm_copy.begin(),sickLm_copy.end(),[&](Sensor* sa,Sensor* sb)->bool
						{
							LineSegment s1_sa(-1,s1,sa);
							LineSegment s1_sb(-1,s1,sb);

							s1_sa.floodStart2EndPathOnItypeBoundary(iType);
							s1_sb.floodStart2EndPathOnItypeBoundary(iType);

							return (s1_sa.itypeBoundaryPath[iType].size() < s1_sb.itypeBoundaryPath[iType].size());
						});

						Sensor* minLenS = *it_min_len;


						int findPosS1_minLens = findItypeEdge(iType,s1,minLenS);
						assert(findPosS1_minLens == -1);
						minLenSensors.insert(minLenS);// the chosen min sensors can not be chosen again.


						LineSegment LH(lineIndex,s1,minLenS);
						LH.floodStart2EndPathOnItypeBoundary(iType);
						iTypeEdgesSet[iType].push_back(LH);
						++ lineIndex;

					}
					else 
					{
						Sensor* s1_nb = *(s1->iTypeSickLandmarks[iType].begin());

						// connect s1 to other sick landmarks except s1_nb
						set<Sensor*>sickLm_copy = sickLm;
						// 
						sickLm_copy.erase(s1_nb);
						sickLm_copy.erase(s1);
						for (auto itMin = minLenSensors.begin(); itMin != minLenSensors.end(); ++ itMin)
						{
							sickLm_copy.erase(*itMin);
						}
						auto it_min_len = min_element(sickLm_copy.begin(),sickLm_copy.end(),[&](Sensor* sa,Sensor* sb)->bool
						{
							LineSegment s1_sa(-1,s1,sa);
							LineSegment s1_sb(-1,s1,sb);

							s1_sa.floodStart2EndPathOnItypeBoundary(iType);
							s1_sb.floodStart2EndPathOnItypeBoundary(iType);

							return (s1_sa.itypeBoundaryPath[iType].size() < s1_sb.itypeBoundaryPath[iType].size());
						});


						Sensor* minLenS = *it_min_len;


						int findPosS1_minLens = findItypeEdge(iType,s1,minLenS);
						assert(findPosS1_minLens == -1);
						minLenSensors.insert(minLenS);// the chosen min sensors can not be chosen again.


						LineSegment LH(lineIndex,s1,minLenS);
						LH.floodStart2EndPathOnItypeBoundary(iType);
						iTypeEdgesSet[iType].push_back(LH);
						++ lineIndex;

						// if the three point form a triangle.
						auto it_find_s1Nb = find(minLenS->iTypeSickLandmarks[iType].begin(),
							minLenS->iTypeSickLandmarks[iType].end(),s1_nb);
						if (it_find_s1Nb != minLenS->iTypeSickLandmarks[iType].end())
						{
							// there should be a triangle.
							int findTriPos = findItypeTriangle(iType,s1,s1_nb,minLenS);
							assert(findTriPos == -1);

							if (findTriPos == -1)
							{
								Triangle newTri(triIndex, s1,s1_nb,minLenS);
								iTypeTriangleSet[iType].push_back(newTri);
								++ triIndex;


								//update edge using times.
								int e12 = findItypeEdge(iType,s1,s1_nb);
								int e23 = findItypeEdge(iType,s1,minLenS);
								int e13 = findItypeEdge(iType,s1_nb,minLenS);

								assert(e12 > -1 && e23 >-1 && e13 >-1);

								iTypeEdgesSet[iType][e12].usedtimes ++;
								iTypeEdgesSet[iType][e23].usedtimes ++;
								iTypeEdgesSet[iType][e13].usedtimes ++;
							}



						}
					}
					

					itypeAttachattachEachSensorWithSickLandmarks(iType);// now every sensor has valid Landmark neighbors.
					itypeAttachEdgesToEachLandmark(iType);
					itypeAttachTriangleToEachLandmark(iType);
					itypeAttachEdgesToTriangles(iType);

					break;
				}
			}

			else
				break;
	
		}
		



		// test if there is a landmark that is no attached valid edges.
		int countSingleLm = 0;
		for (int i=0; i<iType2LandmarkSet[iType].size(); ++i)
		{
			Sensor* s = iType2LandmarkSet[iType][i];
			if (s->iTypeSickLandmarks[iType].size() == 1)
			{
				++countSingleLm;
			}
		}

		cerr<<countSingleLm<<" single landmarks with no valid edges in level "<<iType<<"."<<endl;





		if (countSingleLm > 0 )
		{
			cerr<<"Can not proceed..."<<endl;
			cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<endl;
			exit(EXIT_FAILURE);
		}


		vector<Sensor*>bannedStartEndSensors;
		while (true)
		{
			vector<LineSegment*>sickEdges;
			

			for (int iE = 0; iE < iTypeEdgesSet[iType].size(); ++iE)
			{
				if (iTypeEdgesSet[iType][iE].usedtimes < 2)
				{
					sickEdges.push_back(&iTypeEdgesSet[iType][iE]);
				}
			}

			if (!sickEdges.empty())
			{
				// The 0-edge have the higher priority.
				sort(sickEdges.begin(),sickEdges.end(),[](LineSegment* L1,LineSegment* L2)->bool
				{
					if (L1->usedtimes == L2->usedtimes)
					{
						return (L1->lineImportance > L2->lineImportance);
					}
					else
						return (L1->usedtimes < L2->usedtimes);

				});
				// flooding on the edges to find the holes.
				LineSegment* startEdge = NULL;

				Sensor* srcSensor = NULL;
				Sensor* destSensor = NULL;

				for (auto it_iSickEdge = sickEdges.begin();
					it_iSickEdge != sickEdges.end();++it_iSickEdge)
				{
					startEdge = *it_iSickEdge;
					srcSensor = startEdge->end;
					destSensor = startEdge->start;


					auto itsrc  = find(bannedStartEndSensors.begin(),bannedStartEndSensors.end(),srcSensor);
					auto itdest = find(bannedStartEndSensors.begin(),bannedStartEndSensors.end(),destSensor);

					if (itsrc == bannedStartEndSensors.end() && itdest == bannedStartEndSensors.end())
					{
						break;
					}

				}
				

				// there Edges a sensor attached construct the valid landmark neighbors.
// 				attachEachSensorWithSickLandmarks();// now every sensor has valid Landmark neighbors.
// 				attachEdgesToEachLandmark();
// 				attachTriangleToEachLandmark();
// 				attachEdgesToTriangles();

				vector<Sensor*>holePath;
				//first, remove the relation ship
				srcSensor ->iTypeSickLandmarks[iType].erase(destSensor);
				destSensor->iTypeSickLandmarks[iType].erase(srcSensor);

				srcSensor->floodingOnItypeSickLandmarks(iType,destSensor,holePath);
				// to get a full cycle..
				holePath.insert(holePath.begin(),srcSensor);
				// reattach the landmarks.
				
				itypeAttachattachEachSensorWithSickLandmarks(iType);
				itypeAttachEdgesToEachLandmark(iType);
				itypeAttachTriangleToEachLandmark(iType);
				itypeAttachEdgesToTriangles(iType);

				if (holePath.size() == 3)
				{
					int postri = findItypeTriangle(iType,holePath[0],holePath[1],holePath[2]);
					if (postri >  -1)
					{
						cerr<<"Here is not hole."<<endl;
						bannedStartEndSensors.insert(bannedStartEndSensors.end(),holePath.begin(),holePath.end());
					}
					else
					{
						cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<endl;
						exit(EXIT_FAILURE);
					}
					
				}
				else
				{
					// start filling holes.

					for (int i = 0; i< holePath.size(); ++i)
					{
						Sensor* sStart  =  holePath[i];

						Sensor* s1 = holePath[(i-1+ holePath.size())%holePath.size()];
						Sensor* s2 = holePath[(i+1)%holePath.size()];
						int posS1S2 = findItypeEdge(iType,s1,s2);
						// if the side exist, the line should not be connected, so i give it a very large length.
						if (posS1S2 == -1)
						{
							LineSegment sStart_opposite_Side(-1,s1,s2);

							sStart_opposite_Side.floodStart2EndPathOnItypeBoundary(iType);
							sStart->oppositeSideLength = sStart_opposite_Side.itypeBoundaryPath[iType].size();
						}
						else
							sStart->oppositeSideLength = 0x00FF;

					}

					// find the min length opposite side
					int minLen = 0xFFFF;
					int minPos = -1;
					for (int i = 0; i< holePath.size(); ++i)
					{
						Sensor* s = holePath[i];

						if (s->oppositeSideLength < minLen)
						{
							minLen = s->oppositeSideLength;
							minPos = i;
						}
					}

					// connect the side..
					Sensor* newEdgeStart = holePath[(minPos-1+ holePath.size())%holePath.size()];
					Sensor* newEdgeEnd = holePath[(minPos+1)%holePath.size()];

					int findPos = findItypeEdge(iType,newEdgeStart,newEdgeEnd);

					if (findPos > -1)
					{
						cerr<<"The edge already exist...the hole is wrong..."<<endl;
						cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<endl;
						exit(EXIT_FAILURE);
					}
					else
					{
						// find the two other side, must exist.
						int s2newStart = findItypeEdge(iType,holePath[minPos],newEdgeStart);
						int s2newEnd   = findItypeEdge(iType,holePath[minPos],newEdgeEnd);

						assert(s2newEnd >-1 && s2newStart>-1);

						// there is a new edge...

						LineSegment LH(lineIndex,newEdgeStart,newEdgeEnd);
						LH.floodStart2EndPathOnItypeBoundary(iType);
						iTypeEdgesSet[iType].push_back(LH);
						++ lineIndex;

						//reattachment....
						itypeAttachattachEachSensorWithSickLandmarks(iType);
						itypeAttachEdgesToEachLandmark(iType);
						itypeAttachTriangleToEachLandmark(iType);
						itypeAttachEdgesToTriangles(iType);

						// and there is a new triangle..
						int findTriPos = findItypeTriangle(iType,holePath[minPos],newEdgeStart,newEdgeEnd);
						assert(findTriPos == -1);

						// the 3 edges use times should be modified.
						if (findTriPos == -1)
						{
							//new triangle...
							Triangle newTri(triIndex, holePath[minPos],newEdgeStart,newEdgeEnd);
							iTypeTriangleSet[iType].push_back(newTri);
							++ triIndex;

							//update edge using times.
							int e12 = findItypeEdge(iType,holePath[minPos],newEdgeStart);
							int e23 = findItypeEdge(iType,newEdgeStart,newEdgeEnd);
							int e13 = findItypeEdge(iType,holePath[minPos],newEdgeEnd);

							assert(e12 > -1 && e23 >-1 && e13 >-1);

							iTypeEdgesSet[iType][e12].usedtimes ++;
							iTypeEdgesSet[iType][e23].usedtimes ++;
							iTypeEdgesSet[iType][e13].usedtimes ++;

							assert(iTypeEdgesSet[iType][e12].usedtimes<=2);
							assert(iTypeEdgesSet[iType][e23].usedtimes<=2);
							assert(iTypeEdgesSet[iType][e13].usedtimes<=2);


							if (holePath.size() == 4)
							{
								//there is another triangle...
								Sensor* oppSensor = holePath[(minPos+2)%holePath.size()];
								int findTriPos = findItypeTriangle(iType,oppSensor,newEdgeStart,newEdgeEnd);
								assert(findTriPos == -1);

								if (findTriPos == -1)
								{
									//new triangle...
									Triangle newTri(triIndex, oppSensor,newEdgeStart,newEdgeEnd);
									iTypeTriangleSet[iType].push_back(newTri);
									++ triIndex;

									//update edge using times.
									int e12 = findItypeEdge(iType,oppSensor,newEdgeStart);
									int e23 = findItypeEdge(iType,newEdgeStart,newEdgeEnd);
									int e13 = findItypeEdge(iType,oppSensor,newEdgeEnd);

									assert(e12 > -1 && e23 >-1 && e13 >-1);

									iTypeEdgesSet[iType][e12].usedtimes ++;
									iTypeEdgesSet[iType][e23].usedtimes ++;
									iTypeEdgesSet[iType][e13].usedtimes ++;

									assert(iTypeEdgesSet[iType][e12].usedtimes<=2);
									assert(iTypeEdgesSet[iType][e23].usedtimes<=2);
									assert(iTypeEdgesSet[iType][e13].usedtimes<=2);


								}

							}

							//reattachment....
							itypeAttachattachEachSensorWithSickLandmarks(iType);
							itypeAttachEdgesToEachLandmark(iType);
							itypeAttachTriangleToEachLandmark(iType);
							itypeAttachEdgesToTriangles(iType);

						}



					}


				}
			}
			else 
				break;
		}


		cerr<<iTypeTriangleSet[iType].size()<< " triangles in the ending of the filling hole step in level "<<iType<<"."<<endl;
		cerr<<iTypeEdgesSet[iType].size()<<" edges in the ending of the filling hole step in level "<<iType<<"."<<endl;

		int unAttacchedLM = 0;
		for_each(iType2LandmarkSet[iType].begin(),iType2LandmarkSet[iType].end(),[&](Sensor* s)
		{
			if (s->iTypeAttachedTriangles[iType].empty())
			{
				++ unAttacchedLM;
				cout<<"In level "<<iType<<" , "<<s->index<<" has no attachment."<<endl;
			}
		});

		cout<<unAttacchedLM<<" landmarks in level "<<iType<<endl;

	}

}
void Field::iTypeRicciFlow()
{
	updateTriInfo();//update the triangle neighbors.

	for (int i = 1; i <= MAX_BOUNDARY_TYPE; ++i)
	{
		int iType = i;
		set<Sensor*>iTypeLandmarks;
		for (auto it_iLM =  iType2LandmarkSet[iType].begin();
			      it_iLM != iType2LandmarkSet[iType].end(); ++ it_iLM)
		{
			iTypeLandmarks.insert(*it_iLM);
		}

		// find the minimum id in iTypeLandmarks.
		auto it_miniTypeSensor = max_element(iTypeLandmarks.begin(),iTypeLandmarks.end(),[&](Sensor*s1,Sensor* s2)-> bool
		{
			return (s1->iTypeAttachedTriangles[iType].size() < s2->iTypeAttachedTriangles[iType].size());
		});
		Sensor* miniTypeSensor = *it_miniTypeSensor;

		//Save the removed Surface node.
		iTypeRemovedSurfaceNodes[iType].push_back(miniTypeSensor);
		// remove a node....
		iTypeLandmarks.erase(miniTypeSensor);
		// find the new boundary nodes that are attached with the it_miniTypeSensor
		const int newBoundaryLen = miniTypeSensor->iTypeAttachedTriangles[iType].size();
		set<Sensor*> newBoundarySensors;
		for_each(miniTypeSensor->iTypeAttachedEdges[iType].begin(),miniTypeSensor->iTypeAttachedEdges[iType].end(),[&](LineSegment* L)
		{
			newBoundarySensors.insert(L->start);
			newBoundarySensors.insert(L->end);
		});
		newBoundarySensors.erase(miniTypeSensor);
		assert(newBoundaryLen == newBoundarySensors.size());


		// set the target curvature for every itype landmarks except for miniTypeSensor
		for (auto itS = iTypeLandmarks.begin(); itS !=  iTypeLandmarks.end(); ++ itS)
		{
			Sensor* sLM = *itS;
			auto itFindInBound = find(newBoundarySensors.begin(),newBoundarySensors.end(),sLM);

			if (itFindInBound != newBoundarySensors.end())
			{
				sLM->iTypeTargetCurvature[iType] = 2*PI/newBoundaryLen;
			}
			else
				sLM->iTypeTargetCurvature[iType] = 0.0;

		}
		// set all u...
		for (auto itS = iTypeLandmarks.begin(); itS !=  iTypeLandmarks.end(); ++ itS)
		{
			Sensor* sLM = *itS;
			sLM->iTypeU[iType] = 0.0;
			sLM->iTypeCurrGCurvature[iType] = 0.0;
			sLM->iTypeTheta[iType] = 0.0;
		}
		// find the itype edges...
		set<LineSegment*>iTypeEdges;
		for (int ie = 0; ie < iTypeEdgesSet[iType].size(); ++ ie)
		{
			Sensor* s = iTypeEdgesSet[iType][ie].start;
			Sensor* e = iTypeEdgesSet[iType][ie].end;

			if (s->isBelongToType(iType) && e->isBelongToType(iType) && s != miniTypeSensor && e != miniTypeSensor  )
			{
				iTypeEdges.insert(&iTypeEdgesSet[iType][ie]);
			}
		}
		// find the itype triangles...
		set<Triangle*>iTypeTriangles;
		for (int iT = 0; iT < iTypeTriangleSet[iType].size(); ++ iT)
		{
			Sensor* v1 = iTypeTriangleSet[iType][iT].t1;
			Sensor* v2 = iTypeTriangleSet[iType][iT].t2;
			Sensor* v3 = iTypeTriangleSet[iType][iT].t3;

			if (   v1->isBelongToType(iType) && v2->isBelongToType(iType) && v3->isBelongToType(iType)
				&& v1 != miniTypeSensor && v2 != miniTypeSensor && v3 != miniTypeSensor)
			{
				iTypeTriangles.insert(&iTypeTriangleSet[iType][iT]);
			}
		}
		// start flow...
		while (true)
		{
			// calculate edge length...
			for (auto itE = iTypeEdges.begin(); itE != iTypeEdges.end(); ++ itE)
			{
				LineSegment* LineIJ = *itE;
				Sensor* iSensor = LineIJ->start;
				Sensor* jSensor = LineIJ->end;

				LineIJ->iTypeFlowLength[iType] = exp(iSensor->iTypeU[iType]) + exp(jSensor->iTypeU[iType]);

			}
			// compute the angles...
			for (auto iT = iTypeTriangles.begin(); iT != iTypeTriangles.end(); ++ iT)
			{
				Triangle* Tri_IJK = *iT;
				Sensor* iSensor = Tri_IJK->t1;
				Sensor* jSensor = Tri_IJK->t2;
				Sensor* kSensor = Tri_IJK->t3;

				int edgeIJ = findItypeEdge(iType,iSensor,jSensor);
				int edgeJK = findItypeEdge(iType,jSensor,kSensor);
				int edgeIK = findItypeEdge(iType,iSensor,kSensor);
				// the edges must exist.
				assert(edgeIJ>-1 && edgeJK>-1 && edgeIK>-1);

				long double edgeijflowLen = iTypeEdgesSet[iType][edgeIJ].iTypeFlowLength[iType];
				long double edgejkflowLen = iTypeEdgesSet[iType][edgeJK].iTypeFlowLength[iType];
				long double edgeikflowLen = iTypeEdgesSet[iType][edgeIK].iTypeFlowLength[iType];

				// according to the cosine rule...

				long double lij2 = pow(edgeijflowLen,(long double)2.0);
				long double ljk2 = pow(edgejkflowLen,(long double)2.0);
				long double lik2 = pow(edgeikflowLen,(long double)2.0);

				iSensor->iTypeTheta[iType] += acos((lij2 + lik2 - ljk2)/(2.0* edgeijflowLen * edgeikflowLen ));
				jSensor->iTypeTheta[iType] += acos((lij2 + ljk2 - lik2)/(2.0* edgejkflowLen * edgeijflowLen ));
				kSensor->iTypeTheta[iType] += acos((lik2 + ljk2 - lij2)/(2.0* edgeikflowLen * edgejkflowLen ));

			}

			//update the curvature
			for (auto itS = iTypeLandmarks.begin(); itS !=  iTypeLandmarks.end(); ++ itS)
			{

				Sensor* sLM = *itS;
				if (sLM != miniTypeSensor)
				{
					auto it_findSLM = find(newBoundarySensors.begin(),newBoundarySensors.end(),sLM);
					if (it_findSLM == newBoundarySensors.end())//not belong to the boundary.
					{
						sLM->iTypeCurrGCurvature[iType] = 2* PI - sLM->iTypeTheta[iType];
					}
					else// belong to the boundary.
					{
						sLM->iTypeCurrGCurvature[iType]  =   PI - sLM->iTypeTheta[iType];
					}

				}


			}

			// break the loop if ...
			const long double endEps = 0.00001;
			auto it_max_diffS = max_element(iTypeLandmarks.begin(),iTypeLandmarks.end(),[&](Sensor* s1,Sensor*s2)->bool
			{
				return (fabs(s1->iTypeTargetCurvature[iType] - s1->iTypeCurrGCurvature[iType]) 
					  < fabs(s2->iTypeTargetCurvature[iType] - s2->iTypeCurrGCurvature[iType]));
			});
			double long maxDiff = fabs((*it_max_diffS)->iTypeTargetCurvature[iType] 
									-  (*it_max_diffS)->iTypeCurrGCurvature[iType]);
			cerr<<"Now the maxDiff is "<<maxDiff<<" in "<<iType<<" level."<<endl;
			if (maxDiff < endEps)
			{
				break;
			}

			// update u...
			const long double flowStep = 0.2;
			for_each(iTypeLandmarks.begin(),iTypeLandmarks.end(),[&](Sensor* s)
			{
				s->iTypeU[iType] += flowStep*(s->iTypeTargetCurvature[iType] - s->iTypeCurrGCurvature[iType]);	
				s->iTypeTheta[iType] = 0;// the theta should be reset.
			});


		}

		cerr<<"While loop breaks! "<<endl;

		// Recalculate the angles...
		for (auto iT = iTypeTriangles.begin(); iT != iTypeTriangles.end(); ++ iT)
		{
			Triangle* Tri_IJK = *iT;
			Sensor* iSensor = Tri_IJK->t1;
			Sensor* jSensor = Tri_IJK->t2;
			Sensor* kSensor = Tri_IJK->t3;

			int edgeIJ = findItypeEdge(iType,iSensor,jSensor);
			int edgeJK = findItypeEdge(iType,jSensor,kSensor);
			int edgeIK = findItypeEdge(iType,iSensor,kSensor);
			// the edges must exist.
			assert(edgeIJ>-1 && edgeJK>-1 && edgeIK>-1);


			long double edgeijflowLen = iTypeEdgesSet[iType][edgeIJ].iTypeFlowLength[iType];
			long double edgejkflowLen = iTypeEdgesSet[iType][edgeJK].iTypeFlowLength[iType];
			long double edgeikflowLen = iTypeEdgesSet[iType][edgeIK].iTypeFlowLength[iType];

			// according to the cosine rule...

			long double lij2 = pow(edgeijflowLen,(long double)2.0);
			long double ljk2 = pow(edgejkflowLen,(long double)2.0);
			long double lik2 = pow(edgeikflowLen,(long double)2.0);

			iSensor->iTypeTheta[iType] = acos((lij2 + lik2 - ljk2)/(2.0* edgeijflowLen * edgeikflowLen ));
			jSensor->iTypeTheta[iType] = acos((lij2 + ljk2 - lik2)/(2.0* edgejkflowLen * edgeijflowLen ));
			kSensor->iTypeTheta[iType] = acos((lik2 + ljk2 - lij2)/(2.0* edgeikflowLen * edgejkflowLen ));

			
			Tri_IJK->iTypeV2Theta[iType][iSensor] = iSensor->iTypeTheta[iType];
			Tri_IJK->iTypeV2Theta[iType][jSensor] = jSensor->iTypeTheta[iType];
			Tri_IJK->iTypeV2Theta[iType][kSensor] = kSensor->iTypeTheta[iType];

		}

		// find the boundary edges..
// 		for (int iE = 0; iE < EdgeSet.size(); ++ iE)
// 		{
// 			Sensor* sS = EdgeSet[iE].start;
// 			Sensor* sE = EdgeSet[iE].end;
// 
// 			auto itf_sS = find(newBoundarySensors.begin(),newBoundarySensors.end(),sS);
// 			auto itf_sE = find(newBoundarySensors.begin(),newBoundarySensors.end(),sE);
// 
// 			if (itf_sE != newBoundarySensors.end() && itf_sE != newBoundarySensors.end())
// 			{
// 				EdgeSet[iE].atBoundary = true;
// 			}
// 		}
		// 		for_each(iTypeEdges.begin(),iTypeEdges.end(),[](LineSegment* L)
		// 		{
		// 			cout<<"the Flow length of L: "<<L->flowLength<<endl;
		// 		});
		// 		for_each(iTypeLandmarks.begin(),iTypeLandmarks.end(),[](Sensor* s)
		// 		{
		// 			cout<<" the current curvature of "<<s->index<<" is "<<s->currGCurvature<<endl;
		// 		}); 

		// Planar embedding...
		cout<<"Start planar embedding..."<<endl;
		//clear the state.
		for_each(iTypeTriangles.begin(),iTypeTriangles.end(),[&](Triangle* T)
		{
			T->broadcastCovered = false;
			T->t1->accessed = false;
			T->t2->accessed = false;
			T->t3->accessed = false;
			T->iTypeV2Vmap[iType].clear();
		});

		assert(iTypeTriangles.size() > 1);
		Triangle* seedingTriangle = *(iTypeTriangles.begin());
		Sensor* vs0 = seedingTriangle->t1;
		Sensor* vs1 = seedingTriangle->t2;
		Sensor* vs2 = seedingTriangle->t3;

		// Seeding sequence...
		seedingTriangle->iTypeV2Vmap[iType][vs0] = vs1;
		seedingTriangle->iTypeV2Vmap[iType][vs1] = vs2;
		seedingTriangle->iTypeV2Vmap[iType][vs2] = vs0;

		int t0t1 = findItypeEdge(iType,vs0,vs1);
		int t0t2 = findItypeEdge(iType,vs0,vs2);
		int t1t2 = findItypeEdge(iType,vs1,vs2);
		assert(t0t1 > -1 && t0t2 >-1 && t1t2 > -1);


		long double len_t0t1 = iTypeEdgesSet[iType][t0t1].iTypeFlowLength[iType];
		long double len_t0t2 = iTypeEdgesSet[iType][t0t2].iTypeFlowLength[iType];
		long double t0_theta = seedingTriangle->iTypeV2Theta[iType][vs0];
		// the locations(ricci) of the three points.
		vs0->iTypeRicciLocation[iType] = Point(0,0);
		vs1->iTypeRicciLocation[iType] = Point(len_t0t1,0);
		vs2->iTypeRicciLocation[iType] = Point(len_t0t2*cos(t0_theta),len_t0t2*sin(t0_theta));


		//the first tri is broadcasted.
		vector <Triangle*> broadCastCoverdTris;
		vector <Triangle*> triFrontier;
		triFrontier.push_back(seedingTriangle);

		seedingTriangle->broadcastCovered = true;
		broadCastCoverdTris.push_back(seedingTriangle);

		int accessedTriNumber = 1;

		vs0->accessed = true;
		vs1->accessed = true;
		vs2->accessed = true;

		vector<Sensor*> accessSensorsSet;
		accessSensorsSet.push_back(vs0);
		accessSensorsSet.push_back(vs1);
		accessSensorsSet.push_back(vs2);

		int debugNUM = 0;
		// give every triangle a rotation direction.
		while (accessedTriNumber != iTypeTriangles.size())
		{
// 			cerr<<"Debug "<<debugNUM++<<endl;
// 			if (debugNUM == 0)
// 			{
// 				cerr<<"debug here.."<<endl;
// 			}
			vector<Triangle*>newTriFrontier;
			for (int i=0; i < triFrontier.size(); ++i)
			{
				Triangle* tf = triFrontier[i];
				if (tf->broadcastCovered)
				{
					for (auto   itTriNb = tf->iTypeTriangleNeighbors[iType].begin(); 
						itTriNb != tf->iTypeTriangleNeighbors[iType].end();
						++ itTriNb)
					{
						Triangle* triNb = *itTriNb;
						auto it_findNb = find(iTypeTriangles.begin(),iTypeTriangles.end(),triNb);
						if (! triNb->broadcastCovered && it_findNb != iTypeTriangles.end())
						{
							triNb->broadcastCovered = true;
							broadCastCoverdTris.push_back(triNb);

							newTriFrontier.push_back(triNb);
							++ accessedTriNumber;

							Triangle* frontierNb = triNb->findAItypeNeighborInFrontier(iType, triFrontier);
							LineSegment* sharingEdge = frontierNb->findASharingItypeEdgeWith(iType,triNb);
							Sensor* triNb3_s = triNb->get3rdSensorByEdge(sharingEdge);
							// the edge is oriented in the frontier.
							Sensor* sS = sharingEdge->start;
							Sensor* sE = sharingEdge->end;


							if (frontierNb->iTypeV2Vmap[iType].find(sS)->second == sE)
							{
								triNb->iTypeV2Vmap[iType][sS] = triNb3_s;
								triNb->iTypeV2Vmap[iType][triNb3_s] = sE;
								triNb->iTypeV2Vmap[iType][sE] = sS;
							}
							else
							{
								triNb->iTypeV2Vmap[iType][sE] = triNb3_s;
								triNb->iTypeV2Vmap[iType][triNb3_s] = sS;
								triNb->iTypeV2Vmap[iType][sS] = sE;
							}
						}
					}
				}
			}

			triFrontier.swap(newTriFrontier);
			newTriFrontier.clear();
		}

		triFrontier.clear();
		// recover the state.
		for_each(broadCastCoverdTris.begin(),broadCastCoverdTris.end(),[](Triangle* T1)
		{
			T1->broadcastCovered = false;
		});

		// give every vertex its ricci location.
		debugNUM = 0;
		int accessedSensorNum = 3;
		while (accessedSensorNum != iTypeLandmarks.size())
		{
			for (auto it_iLms = iTypeLandmarks.begin(); it_iLms != iTypeLandmarks.end(); ++ it_iLms)
			{
				// 				cerr<<"Debug num..."<<debugNUM++<<endl;
				// 				if (debugNUM == 622)
				// 				{
				// 					cerr<<"debug here..."<<endl;
				// 				}
				Sensor* iLms = *it_iLms;
				if (! iLms->accessed)
				{
					LineSegment* accessedEdge = iLms->findAnItypeAccessedEdgeInAttachedTriangle(iType);
					if (accessedEdge != NULL)
					{
						Triangle* T = findItypeTrianglePointer(iType,accessedEdge->start,accessedEdge->end,iLms);
						//cout<<"Distance diff"<<accessedEdge->start->ricciLocation.distance(accessedEdge->end->ricciLocation)
						// - accessedEdge->flowLength<<endl;

						LineSegment* start_ilms = findItypeEdgePointer(iType,accessedEdge->start,iLms);
						LineSegment* end_ilms   = findItypeEdgePointer(iType,accessedEdge->end, iLms);

						long double len_start_ilms = start_ilms->iTypeFlowLength[iType];
						long double len_end_ilms   = end_ilms->iTypeFlowLength[iType];

						if (T->iTypeV2Vmap[iType][accessedEdge->start] == accessedEdge->end)
						{
							iLms->calPosWithCircleIntersecCircleItype(iType,
								accessedEdge->end->iTypeRicciLocation[iType],len_end_ilms,
								accessedEdge->start->iTypeRicciLocation[iType],len_start_ilms);

							iLms->accessed = true;
							accessSensorsSet.push_back(iLms);

							++ accessedSensorNum;
						}
						else
						{
							iLms->calPosWithCircleIntersecCircleItype(iType,
								accessedEdge->start->iTypeRicciLocation[iType],len_start_ilms,
								accessedEdge->end->iTypeRicciLocation[iType],len_end_ilms);

							iLms->accessed = true;
							accessSensorsSet.push_back(iLms);

							++ accessedSensorNum;
						}
					}
				}
			}
		}

		for_each(accessSensorsSet.begin(),accessSensorsSet.end(),[](Sensor* s)
		{
			s->accessed = false;
		});

		cout<<"End of planarization..."<<endl;
	}

}
void Field::ricciFlow()
{

	for (int i = 0; i <= MAX_BOUNDARY_TYPE; ++i)
	{
 		set<Sensor*>iTypeLandmarks;
		for (int ilm = 0;ilm < LandmarkSet.size(); ++ilm)
		{
			Sensor* itypelm = LandmarkSet[ilm];
			if (itypelm->boundarySensorType == i)
			{
				iTypeLandmarks.insert(itypelm);
			}
		}
	
	
		// find the minimum id in iTypeLandmarks.
		auto it_miniTypeSensor = max_element(iTypeLandmarks.begin(),iTypeLandmarks.end(),[](Sensor*s1,Sensor* s2)-> bool
		{
			return (s1->attachedTriangles.size() < s2->attachedTriangles.size());
		});
		Sensor* miniTypeSensor = *it_miniTypeSensor;

		//Save the removed Surface node.
		removedSurfaceNodes.push_back(miniTypeSensor);
		// remove a node....
		iTypeLandmarks.erase(miniTypeSensor);
		// find the new boundary nodes that are attached with the it_miniTypeSensor
		const int newBoundaryLen = miniTypeSensor->attachedTriangles.size();
		set<Sensor*> newBoundarySensors;
		for_each(miniTypeSensor->attachedEdges.begin(),miniTypeSensor->attachedEdges.end(),[&](LineSegment* L)
		{
			newBoundarySensors.insert(L->start);
			newBoundarySensors.insert(L->end);
		});
		newBoundarySensors.erase(miniTypeSensor);
		assert(newBoundaryLen == newBoundarySensors.size());

		
		// set the target curvature for every itype landmarks except for miniTypeSensor
		for (auto itS = iTypeLandmarks.begin(); itS !=  iTypeLandmarks.end(); ++ itS)
		{
			Sensor* sLM = *itS;
			auto itFindInBound = find(newBoundarySensors.begin(),newBoundarySensors.end(),sLM);

			if (itFindInBound != newBoundarySensors.end())
			{
				sLM->targetCurvature = 2*PI/newBoundaryLen;
			}
			else
				sLM->targetCurvature = 0.0;
			
		}
		// set all u...
		for (auto itS = iTypeLandmarks.begin(); itS !=  iTypeLandmarks.end(); ++ itS)
		{
			Sensor* sLM = *itS;
			sLM->u = 0.0;
			sLM->currGCurvature = 0.0;
			sLM->theta = 0.0;
		}
		// find the itype edges...
		set<LineSegment*>iTypeEdges;
		for (int ie = 0; ie < EdgeSet.size(); ++ ie)
		{
			Sensor* s = EdgeSet[ie].start;
			Sensor* e = EdgeSet[ie].end;

			if (s->boundarySensorType == i && s != miniTypeSensor && e != miniTypeSensor  )
			{
				iTypeEdges.insert(&EdgeSet[ie]);
			}
		}
		// find the itype triangles...
		set<Triangle*>iTypeTriangles;
		for (int iT = 0; iT < TriangleSet.size(); ++ iT)
		{
			Sensor* v1 = TriangleSet[iT].t1;
			Sensor* v2 = TriangleSet[iT].t2;
			Sensor* v3 = TriangleSet[iT].t3;

			if (v1->boundarySensorType == i && v1 != miniTypeSensor && v2 != miniTypeSensor && v3 != miniTypeSensor)
			{
				iTypeTriangles.insert(&TriangleSet[iT]);
			}
		}
		// start flow...
		while (true)
		{
			// calculate edge length...
			for (auto itE = iTypeEdges.begin(); itE != iTypeEdges.end(); ++ itE)
			{
				LineSegment* LineIJ = *itE;
				Sensor* iSensor = LineIJ->start;
				Sensor* jSensor = LineIJ->end;

				LineIJ->flowLength = exp(iSensor->u) + exp(jSensor->u);

			}
			// compute the angles...
			for (auto iT = iTypeTriangles.begin(); iT != iTypeTriangles.end(); ++ iT)
			{
				Triangle* Tri_IJK = *iT;
				Sensor* iSensor = Tri_IJK->t1;
				Sensor* jSensor = Tri_IJK->t2;
				Sensor* kSensor = Tri_IJK->t3;

				int edgeIJ = findEdge(iSensor,jSensor);
				int edgeJK = findEdge(jSensor,kSensor);
				int edgeIK = findEdge(iSensor,kSensor);
				// the edges must exist.
				assert(edgeIJ>-1 && edgeJK>-1 && edgeIK>-1);
				
				long double edgeijflowLen = EdgeSet[edgeIJ].flowLength;
				long double edgejkflowLen = EdgeSet[edgeJK].flowLength;
				long double edgeikflowLen = EdgeSet[edgeIK].flowLength;

				// according to the cosine rule...

				long double lij2 = pow(edgeijflowLen,(long double)2.0);
				long double ljk2 = pow(edgejkflowLen,(long double)2.0);
				long double lik2 = pow(edgeikflowLen,(long double)2.0);

				iSensor->theta += acos((lij2 + lik2 - ljk2)/(2.0* edgeijflowLen * edgeikflowLen ));
				jSensor->theta += acos((lij2 + ljk2 - lik2)/(2.0* edgejkflowLen * edgeijflowLen ));
				kSensor->theta += acos((lik2 + ljk2 - lij2)/(2.0* edgeikflowLen * edgejkflowLen ));

			}

			//update the curvature
			for (auto itS = iTypeLandmarks.begin(); itS !=  iTypeLandmarks.end(); ++ itS)
			{
			
				Sensor* sLM = *itS;
				if (sLM != miniTypeSensor)
				{
					auto it_findSLM = find(newBoundarySensors.begin(),newBoundarySensors.end(),sLM);
					if (it_findSLM == newBoundarySensors.end())//not belong to the boundary.
					{
						sLM->currGCurvature = 2* PI - sLM->theta;
					}
					else// belong to the boundary.
					{
						sLM->currGCurvature = PI - sLM->theta;
					}
					
				}

				
			}

			// break the loop if ...
			const long double endEps = 0.001;
			auto it_max_diffS = max_element(iTypeLandmarks.begin(),iTypeLandmarks.end(),[&](Sensor* s1,Sensor*s2)->bool
			{
				return (fabs(s1->targetCurvature - s1->currGCurvature) < fabs(s2->targetCurvature - s2->currGCurvature));
			});
			double long maxDiff = fabs((*it_max_diffS)->targetCurvature - (*it_max_diffS)->currGCurvature);
			cerr<<"Now the maxDiff is "<<maxDiff<<endl;
			if (maxDiff < endEps)
			{
				break;
			}

			// update u...
			const double flowStep = 0.2;
			for_each(iTypeLandmarks.begin(),iTypeLandmarks.end(),[&](Sensor* s)
			{
				s->u += flowStep*(s->targetCurvature - s->currGCurvature);	
				s->theta = 0;// the theta should be reset.
			});
			
			
		}

		cerr<<"While loop breaks! "<<endl;

		// Recalculate the angles...
		for (auto iT = iTypeTriangles.begin(); iT != iTypeTriangles.end(); ++ iT)
		{
			Triangle* Tri_IJK = *iT;
			Sensor* iSensor = Tri_IJK->t1;
			Sensor* jSensor = Tri_IJK->t2;
			Sensor* kSensor = Tri_IJK->t3;

			int edgeIJ = findEdge(iSensor,jSensor);
			int edgeJK = findEdge(jSensor,kSensor);
			int edgeIK = findEdge(iSensor,kSensor);
			// the edges must exist.
			assert(edgeIJ>-1 && edgeJK>-1 && edgeIK>-1);

			long double edgeijflowLen = EdgeSet[edgeIJ].flowLength;
			long double edgejkflowLen = EdgeSet[edgeJK].flowLength;
			long double edgeikflowLen = EdgeSet[edgeIK].flowLength;

			// according to the cosine rule...

			long double lij2 = pow(edgeijflowLen,(long double)2.0);
			long double ljk2 = pow(edgejkflowLen,(long double)2.0);
			long double lik2 = pow(edgeikflowLen,(long double)2.0);

			iSensor->theta = acos((lij2 + lik2 - ljk2)/(2.0* edgeijflowLen * edgeikflowLen ));
			jSensor->theta = acos((lij2 + ljk2 - lik2)/(2.0* edgejkflowLen * edgeijflowLen ));
			kSensor->theta = acos((lik2 + ljk2 - lij2)/(2.0* edgeikflowLen * edgejkflowLen ));

			Tri_IJK->v2theta[iSensor] = iSensor->theta;
			Tri_IJK->v2theta[jSensor] = jSensor->theta;
			Tri_IJK->v2theta[kSensor] = kSensor->theta;

		}

		// find the boundary edges..
		for (int iE = 0; iE < EdgeSet.size(); ++ iE)
		{
			Sensor* sS = EdgeSet[iE].start;
			Sensor* sE = EdgeSet[iE].end;

			auto itf_sS = find(newBoundarySensors.begin(),newBoundarySensors.end(),sS);
			auto itf_sE = find(newBoundarySensors.begin(),newBoundarySensors.end(),sE);

			if (itf_sE != newBoundarySensors.end() && itf_sE != newBoundarySensors.end())
			{
				EdgeSet[iE].atBoundary = true;
			}
		}
// 		for_each(iTypeEdges.begin(),iTypeEdges.end(),[](LineSegment* L)
// 		{
// 			cout<<"the Flow length of L: "<<L->flowLength<<endl;
// 		});
// 		for_each(iTypeLandmarks.begin(),iTypeLandmarks.end(),[](Sensor* s)
// 		{
// 			cout<<" the current curvature of "<<s->index<<" is "<<s->currGCurvature<<endl;
// 		}); 

		// Planar embedding...
		cout<<"Start planar embedding..."<<endl;
		//clear the state.
		for_each(iTypeTriangles.begin(),iTypeTriangles.end(),[](Triangle* T)
		{
			T->broadcastCovered = false;
			T->t1->accessed = false;
			T->t2->accessed = false;
			T->t3->accessed = false;
			T->v2vMap.clear();
		});

		assert(iTypeTriangles.size() > 1);
		Triangle* seedingTriangle = *(iTypeTriangles.begin());
		Sensor* vs0 = seedingTriangle->t1;
		Sensor* vs1 = seedingTriangle->t2;
		Sensor* vs2 = seedingTriangle->t3;

		// Seeding sequence...
		seedingTriangle->v2vMap[vs0] = vs1;
		seedingTriangle->v2vMap[vs1] = vs2;
		seedingTriangle->v2vMap[vs2] = vs0;

		int t0t1 = findEdge(vs0,vs1);
		int t0t2 = findEdge(vs0,vs2);
		int t1t2 = findEdge(vs1,vs2);
		assert(t0t1 > -1 && t0t2 >-1 && t1t2 > -1);


		long double len_t0t1 = EdgeSet[t0t1].flowLength;
		long double len_t0t2 = EdgeSet[t0t2].flowLength;
		long double t0_theta = seedingTriangle->v2theta[vs0];
		// the locations(ricci) of the three points.
		vs0->ricciLocation = Point(0,0);
		vs1->ricciLocation = Point(len_t0t1,0);
		vs2->ricciLocation = Point(len_t0t2*cos(t0_theta),len_t0t2*sin(t0_theta));

		
		//the first tri is broadcasted.
		vector<Triangle*> triFrontier;
		triFrontier.push_back(seedingTriangle);
		seedingTriangle->broadcastCovered = true;
		int accessedTriNumber = 1;

		vs0->accessed = true;
		vs1->accessed = true;
		vs2->accessed = true;

		int debugNUM = 0;
		// give every triangle a rotation direction.
		while (accessedTriNumber != iTypeTriangles.size())
		{
			cerr<<"Debug "<<debugNUM++<<endl;
			if (debugNUM == 0)
			{
				cerr<<"debug here.."<<endl;
			}
			vector<Triangle*>newTriFrontier;
			for (int i=0; i < triFrontier.size(); ++i)
			{
				Triangle* tf = triFrontier[i];
				if (tf->broadcastCovered)
				{
					for (auto   itTriNb = tf->TriangleNeighbors.begin(); 
						itTriNb != tf->TriangleNeighbors.end();
						++ itTriNb)
					{
						Triangle* triNb = *itTriNb;
						auto it_findNb = find(iTypeTriangles.begin(),iTypeTriangles.end(),triNb);
						if (! triNb->broadcastCovered && it_findNb != iTypeTriangles.end())
						{
							triNb->broadcastCovered = true;
							newTriFrontier.push_back(triNb);
							++ accessedTriNumber;

							Triangle* frontierNb = triNb->findANeighborInFrontier(triFrontier);
							LineSegment* sharingEdge = frontierNb->findASharingEdgeWith(triNb);
							Sensor* triNb3_s = triNb->get3rdSensorByEdge(sharingEdge);
							// the edge is oriented in the frontier.
							Sensor* sS = sharingEdge->start;
							Sensor* sE = sharingEdge->end;


							if (frontierNb->v2vMap.find(sS)->second == sE)
							{
								triNb->v2vMap[sS] = triNb3_s;
								triNb->v2vMap[triNb3_s] = sE;
								triNb->v2vMap[sE] = sS;
							}
							else
							{
								triNb->v2vMap[sE] = triNb3_s;
								triNb->v2vMap[triNb3_s] = sS;
								triNb->v2vMap[sS] = sE;
							}
						}
					}
				}
			}

			triFrontier.swap(newTriFrontier);
			newTriFrontier.clear();
		}

		triFrontier.clear();


		// give every vertex its ricci location.
		debugNUM = 0;
		int accessedSensorNum = 3;
		while (accessedSensorNum != iTypeLandmarks.size())
		{
			
			for (auto it_iLms = iTypeLandmarks.begin(); it_iLms != iTypeLandmarks.end(); ++ it_iLms)
			{
// 				cerr<<"Debug num..."<<debugNUM++<<endl;
// 				if (debugNUM == 622)
// 				{
// 					cerr<<"debug here..."<<endl;
// 				}
				Sensor* iLms = *it_iLms;
				if (! iLms->accessed)
				{
					LineSegment* accessedEdge = iLms->findAnAccessedEdgeInAttachedTriangle();
					if (accessedEdge != NULL)
					{
						Triangle* T = findTrianglePointer(accessedEdge->start,accessedEdge->end,iLms);
						//cout<<"Distance diff"<<accessedEdge->start->ricciLocation.distance(accessedEdge->end->ricciLocation)
						// - accessedEdge->flowLength<<endl;

						LineSegment* start_ilms = findEdgePointer(accessedEdge->start,iLms);
						LineSegment* end_ilms   = findEdgePointer(accessedEdge->end, iLms);

						long double len_start_ilms = start_ilms->flowLength;
						long double len_end_ilms   = end_ilms->flowLength;

						if (T->v2vMap[accessedEdge->start] == accessedEdge->end)
						{
							iLms->calPosWithCircleIntersecCircle(accessedEdge->end->ricciLocation,len_end_ilms,
								accessedEdge->start->ricciLocation,len_start_ilms);
							
							iLms->accessed = true;
							++ accessedSensorNum;
						}
						else
						{
							iLms->calPosWithCircleIntersecCircle(accessedEdge->start->ricciLocation,len_start_ilms,
								accessedEdge->end->ricciLocation,len_end_ilms);

							iLms->accessed = true;
							++ accessedSensorNum;
						}
					}
				}
			}
		}
		


		cout<<"End of planarization..."<<endl;
	}
	


}

void Field::stereographicProjection()
{
	//first, map the ricci locations to a sphere surface.
// 	for (int ilm = 0; ilm < LandmarkSet.size(); ++ ilm)
// 	{
// 		Sensor* lm = LandmarkSet[ilm];
// 		auto itRemoved = find(removedSurfaceNodes.begin(),removedSurfaceNodes.end(),lm);
// 		if (itRemoved != removedSurfaceNodes.end())
// 		{
// 			Sensor* sRemoved = *itRemoved;
// 			sRemoved->sphereLocation = Point(0,0,1);
// 		}
// 		else
// 			lm->ricci2Sphere();
// 	}
	// map the sphere to a unit disk.
// 	for (int ilm = 0; ilm < LandmarkSet.size(); ++ ilm)
// 	{
// 		Sensor* lm = LandmarkSet[ilm];
// 		auto itRemoved = find(removedSurfaceNodes.begin(),removedSurfaceNodes.end(),lm);
// 		if (itRemoved != removedSurfaceNodes.end())
// 		{
// 			Sensor* sRemoved = *itRemoved;
// 			sRemoved->unitDiskLocation = Point(100,100);
// 		}
// 		else
// 			lm->sphere2Unitdisk();
// 	}

	// mobius transformation.
	//Clear the broadcast state.
	for_each(LandmarkSet.begin(),LandmarkSet.end(),[](Sensor* s)
	{
		s->broadcastCovered = 0;
	});
	//build edge neighbors...
	buildEdgeAttachedLandmarkNb();
	//flood on the landmark nbs to get the longest shortest distance landmark.
	for (int i = 0; i <= MAX_BOUNDARY_TYPE; ++i)
	{
		set<Sensor*>iTypeLandmarks;
		for (int ilm = 0;ilm < LandmarkSet.size(); ++ilm)
		{
			Sensor* itypelm = LandmarkSet[ilm];
			if (itypelm->boundarySensorType == i)
			{
				iTypeLandmarks.insert(itypelm);
			}
		}
		Sensor* iRemoved = removedSurfaceNodes[i];
		iTypeLandmarks.erase(iRemoved);

		Sensor* maxLenSensor = NULL;
		int maxLen = 0;
		int flen = 0;
		Sensor* v3 = *(iRemoved->edgeAttachedLandmarkNb.begin());
		
		for_each(iTypeLandmarks.begin(),iTypeLandmarks.end(),[&](Sensor* sLm)
		{
			flen = v3->floodOnEdgeAttachedLandmarkNbs(sLm);
			if (flen > maxLen)
			{
				maxLenSensor = sLm;
			}
		});
		// find the path..
		assert(maxLenSensor != v3 && maxLenSensor != iRemoved);
		vector<Sensor*> longestPath;
		v3->floodOnEdgeAttachedLandmarkNbsRoutes(maxLenSensor,longestPath);
		// v1 longestPath.back, v2 longestPath.middle, v3 iRemoved.
		assert(longestPath.size() > 2);
		// perform conformal mapping
		Sensor* v1 = longestPath.back();
		Sensor* v2 = longestPath[(int)(longestPath.size()-1)/2];

		complex<long double>z1(v1->ricciLocation.x,v1->ricciLocation.y);
		complex<long double>z2(v2->ricciLocation.x,v2->ricciLocation.y);
		complex<long double>z3(v3->ricciLocation.x,v3->ricciLocation.y);

		for (auto itS = iTypeLandmarks.begin(); itS != iTypeLandmarks.end(); ++ itS)
		{
			Sensor* s = *itS;
			if (s!=v1 && s!= v2 && s!= v3)
			{
				complex<long double>mappedLocation;
				complex<long double>originalLocation(s->ricciLocation.x,s->ricciLocation.y);
				mappedLocation = (originalLocation - z1) * (z2 - z3)/((originalLocation - z3)*(z2 - z1));
				s->conformalLocation.x = mappedLocation.real();
				s->conformalLocation.y = mappedLocation.imag();
			}
			
		}
		// For v1, v2 and v3
		long double nan = std::numeric_limits<double>::quiet_NaN();
		v1->conformalLocation = Point(0,0);
		v2->conformalLocation = Point(1,0);
		v3->conformalLocation = Point(nan,nan);
		// sphere location
		for (auto itS = iTypeLandmarks.begin(); itS != iTypeLandmarks.end(); ++ itS)
		{
			Sensor* s = *itS;
			if ( s!= v3)
			{
				s->conformal2Sphere();
				
			}

		}
		v3->sphereLocation = Point(0,0,1);
		iRemoved->sphereLocation = Point(0,0,1);// should be ...
	}

}
void Sensor::floodOnEdgeAttachedLandmarkNbsRoutes(Sensor* destSensor, vector<Sensor*>& hopPath)
{
	// destSensor should not be null.
	assert(destSensor);
	// This node should be in the same boundary with destSensor.
	assert(boundarySensorType == destSensor->boundarySensorType);
	// the hop path should be empty.
	assert(hopPath.empty());

	int thisBoundaryType = boundarySensorType;
	int i;
	int hopCount = 0;
	set<int> touchedNodes;

	//trivial case
	if (destSensor == this)
		return;

	//the path is in the cache
	if (spMapOnEdgeAttachedLandmarks.find(destSensor->index) != spMapOnEdgeAttachedLandmarks.end())
	{
		Sensor* curNode = this;

		while (curNode != destSensor)
		{
			Sensor *curParent = &globalField->sensorPool[curNode->spMapOnEdgeAttachedLandmarks[destSensor->index].second];
			hopPath.push_back(curParent);
			curNode = curParent;
		};

		return;
	}

	vector<Sensor*> *frontier = new vector<Sensor*>;
	frontier->push_back(this);
	touchedNodes.insert(index);
	broadcastCovered = 1;

	while (destSensor->broadcastCovered != 1)
	{
		vector<Sensor*> *newFrontier = new vector<Sensor*>;

		for(i = 0; i < (int)(frontier->size()); i++)
		{
			Sensor *s = (*frontier)[i];

			if (s->broadcastCovered == 1)
			{
				for (auto it_sNb = s->edgeAttachedLandmarkNb.begin(); it_sNb != s->edgeAttachedLandmarkNb.end(); ++ it_sNb)
				{
					Sensor *n = *it_sNb;
					// this condition limited the flooding on the boundary.
					if (n->broadcastCovered == 0 && thisBoundaryType == n->boundarySensorType)
					{
						n->spMapOnEdgeAttachedLandmarks.insert(DEST2SP::value_type(index, SPLEN_NEXTHOP(hopCount + 1, s->index)));
						n->broadcastCovered = 1;
						touchedNodes.insert(n->index);
						newFrontier->push_back(n);
					}
				}
			}
		}

		delete frontier;
		frontier = newFrontier;

		hopCount++;
	}

	delete frontier;

	Sensor* curNode = destSensor;

	do
	{
		hopPath.insert(hopPath.begin(), curNode);
		Sensor *curParent = &globalField->sensorPool[curNode->spMapOnEdgeAttachedLandmarks[index].second];
		//also remember the reverse route
		curParent->spMapOnEdgeAttachedLandmarks.insert(DEST2SP::value_type(destSensor->index, SPLEN_NEXTHOP(hopPath.size(), curNode->index)));
		curNode = curParent;
	}
	while (curNode != this);

	for(set<int>::const_iterator iter = touchedNodes.begin(); iter != touchedNodes.end(); ++iter)
	{
		globalField->sensorPool[*iter].broadcastCovered = 0;
		globalField->sensorPool[*iter].parent = NULL;
	}


}
int Sensor::floodOnEdgeAttachedLandmarkNbs(Sensor* destSensor)
{
	
  	// destSensor should not be null.
	assert(destSensor);
	// This node should be in the same boundary with destSensor.
	assert(boundarySensorType == destSensor->boundarySensorType);

	int thisBoundaryType = boundarySensorType;
	int i;
	int hopCount = 0;
	set<int> touchedNodes;

	int pathLength = 0;

	//trivial case
	if (destSensor == this)
		return 0;

	//the path is in the cache
	if (spMapOnEdgeAttachedLandmarks.find(destSensor->index) != spMapOnEdgeAttachedLandmarks.end())
	{
		return spMapOnEdgeAttachedLandmarks[destSensor->index].first;
	}

	vector<Sensor*> *frontier = new vector<Sensor*>;
	frontier->push_back(this);
	touchedNodes.insert(index);
	broadcastCovered = 1;

	while (destSensor->broadcastCovered != 1)
	{
		vector<Sensor*> *newFrontier = new vector<Sensor*>;

		for(i = 0; i < (int)(frontier->size()); i++)
		{
			Sensor *s = (*frontier)[i];

			if (s->broadcastCovered == 1)
			{
				for (auto it_sNb = s->edgeAttachedLandmarkNb.begin(); it_sNb != s->edgeAttachedLandmarkNb.end(); ++ it_sNb)
				{
					Sensor *n = *it_sNb;
					// this condition limited the flooding on the boundary.
					if (n->broadcastCovered == 0 && thisBoundaryType == n->boundarySensorType)
					{
						n->spMapOnEdgeAttachedLandmarks.insert(DEST2SP::value_type(index, SPLEN_NEXTHOP(hopCount + 1, s->index)));
						n->broadcastCovered = 1;
						touchedNodes.insert(n->index);
						newFrontier->push_back(n);
					}
				}
			}
		}

		delete frontier;
		frontier = newFrontier;

		hopCount++;
	}

	delete frontier;

	Sensor* curNode = destSensor;

	do
	{
		// boundPath.insert(boundPath.begin(), curNode);
		++ pathLength;
		Sensor *curParent = &globalField->sensorPool[curNode->spMapOnEdgeAttachedLandmarks[index].second];
		//also remember the reverse route
		curParent->spMapOnEdgeAttachedLandmarks.insert(DEST2SP::value_type(destSensor->index, SPLEN_NEXTHOP(pathLength, curNode->index)));
		curNode = curParent;
	}
	while (curNode != this);

	for(set<int>::const_iterator iter = touchedNodes.begin(); iter != touchedNodes.end(); ++iter)
	{
		globalField->sensorPool[*iter].broadcastCovered = 0;
		globalField->sensorPool[*iter].parent = NULL;
	}


	return pathLength;

}
void Sensor::sphere2Unitdisk()
{
// 	if (sphereLocation.z == 1.0)
// 	{
// 		cerr<<"Mapping from sphere to unit disk failed..."<<endl;
// 		cerr<<"The divider would be 0!---fatal error!"<<endl;
// 		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
// 		exit(EXIT_FAILURE);
// 	}
// 	
// 	unitDiskLocation.x = sphereLocation.x / (1-sphereLocation.z);
// 	unitDiskLocation.y = sphereLocation.y / (1-sphereLocation.z);

}
void Sensor::conformal2Sphere()
{
	long double denom =  1.0 + pow(conformalLocation.x,(long double)2.0) + pow(conformalLocation.y,(long double)2.0);
	long double znume = -1.0 + pow(conformalLocation.x,(long double)2.0) + pow(conformalLocation.y,(long double)2.0);

	sphereLocation.x = (2 * conformalLocation.x) / denom;
	sphereLocation.y = (2 * conformalLocation.y) / denom;
	sphereLocation.z = znume / denom;
}
void Field::dumpSphereLocation(string file_out_name)
{
	ofstream file_out;
	file_out.open(file_out_name.c_str(),ofstream::trunc);
	if ( !file_out )
	{
		cerr << "\n";
		cerr << "  Cannot open the output file \"" << file_out_name << "\".\n";
		cerr << "Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<"\n";
		exit(EXIT_FAILURE);
	}
	for (int i = 0; i < nSensors; ++ i)
	{
		Sensor* s = &sensorPool[i];
		file_out<< s->index<<" "<<s->landmark<<" ";
		file_out<<s->sphereLocation.x<<" "<<s->sphereLocation.y<<" "<<s->sphereLocation.z<<" ";
		file_out<<s->boundarySensorType<<endl;
	}

	file_out.close();
}
void Field::dumpLayersOfTriangle(string file_out_name)
{
	ofstream file_out;
	file_out.open(file_out_name.c_str(),ofstream::trunc);
	if ( !file_out )
	{
		cerr << "\n";
		cerr << "  Cannot open the output file \"" << file_out_name << "\".\n";
		cerr << "Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<"\n";
		exit(EXIT_FAILURE);
	}
	for (int iType = 1; iType < MAX_BOUNDARY_TYPE; ++ iType)
	{
		for (auto i_triangle = iTypeTriangleSet[iType].begin(); i_triangle != iTypeTriangleSet[iType].end();
			++ i_triangle)
		{
			Sensor* s1 = (*i_triangle).t1;
			Sensor* s2 = (*i_triangle).t2;
			Sensor* s3 = (*i_triangle).t3;
			
			file_out<<iType<<" ";
			file_out<<s1->index<<" ";
			file_out<<s2->index<<" ";
			file_out<<s3->index;
			file_out<<endl;
		}
// 		Sensor* s = &sensorPool[i];
// 		file_out << s->index<<" ";
// 		file_out<<s->landmark<<" "<<s->ricciLocation.x<<" "<<s->ricciLocation.y<<" "<<s->boundarySensorType<<endl;
	}

	file_out.close();
}
void Field::dumpRicciLocation(string file_out_name)
{
	ofstream file_out;
	file_out.open(file_out_name.c_str(),ofstream::trunc);
	if ( !file_out )
	{
		cerr << "\n";
		cerr << "  Cannot open the output file \"" << file_out_name << "\".\n";
		cerr << "Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<"\n";
		exit(EXIT_FAILURE);
	}
	for (int i = 0; i < nSensors; ++ i)
	{
		Sensor* s = &sensorPool[i];
		file_out << s->index<<" "<<s->landmark<<" "<<s->ricciLocation.x<<" "<<s->ricciLocation.y<<" "<<s->boundarySensorType<<endl;
	}

	file_out.close();

}
int Field::test3ItypeEdges1usedtimeHole(const int itype)
{
	int n3Edges1usedtimeHole = 0;
	int removedTrisCounter = 0;
	for (int iE = 0; iE < iTypeEdgesSet[itype].size(); ++ iE)
	{
		for ( int jE = iE + 1; jE < iTypeEdgesSet[itype].size(); ++ jE )
		{
			if (iTypeEdgesSet[itype][iE].TheSameLine(&(iTypeEdgesSet[itype][jE])))
			{
				// the two line should not be the same line
				cerr<<"Same lines in the EdgeSet."<<endl;
				cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<endl;
				exit(EXIT_FAILURE);
			}
			else
			{
				//the two edges should sharing a vertex and the usedTimes have special relationship.
				Sensor* sharedV = iTypeEdgesSet[itype][iE].sharedVertice(iTypeEdgesSet[itype][jE]);
				if (sharedV != NULL && iTypeEdgesSet[itype][iE].usedtimes == 1 && iTypeEdgesSet[itype][jE].usedtimes == 1)
				{
					vector<Sensor*>assumedTriVertex;
					assumedTriVertex.push_back(iTypeEdgesSet[itype][iE].start);
					assumedTriVertex.push_back(iTypeEdgesSet[itype][iE].end);

					assumedTriVertex.push_back(iTypeEdgesSet[itype][jE].start);
					assumedTriVertex.push_back(iTypeEdgesSet[itype][jE].end);


					assumedTriVertex.erase(remove(assumedTriVertex.begin(),assumedTriVertex.end(),sharedV),assumedTriVertex.end());
					assert(assumedTriVertex.size() == 2);// size should be 2

					int linPos = findItypeEdge(itype,assumedTriVertex[0],assumedTriVertex[1]);
					if (linPos > -1 && iTypeEdgesSet[itype][linPos].usedtimes == 1)
					{
						assert(iTypeEdgesSet[itype][linPos].itypeEdgeAttachedTriangles[itype].size() == 1);
						int triPos = findItypeTriangle(itype,sharedV,assumedTriVertex[0],assumedTriVertex[1]);
						if (triPos == -1)
						{
							++ n3Edges1usedtimeHole;
						}
						Triangle TTT(-1,sharedV,assumedTriVertex[0],assumedTriVertex[1]);
						auto it_removed = find_if(iTypeRemovedTris[itype].begin(),iTypeRemovedTris[itype].end(),[&](Triangle& tInS)
						{
							return (tInS == TTT);
						});
						if (it_removed!=iTypeRemovedTris[itype].end())
						{
							cerr<<"Yes...";
						}

					}
				}
			}
		}
	}

	return n3Edges1usedtimeHole;
}
int Field::test3Edges1usedtimeHole()
{
	int n3Edges1usedtimeHole = 0;
	int removedTrisCounter = 0;
	for (int iE = 0; iE < EdgeSet.size(); ++ iE)
	{
		for ( int jE = iE + 1; jE < EdgeSet.size(); ++ jE )
		{
			if (EdgeSet[iE].TheSameLine(&(EdgeSet[jE])))
			{
				// the two line should not be the same line
				cerr<<"Same lines in the EdgeSet."<<endl;
				cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<endl;
				exit(EXIT_FAILURE);
			}
			else
			{
				//the two edges should sharing a vertex and the usedTimes have special relationship.
				Sensor* sharedV = EdgeSet[iE].sharedVertice(EdgeSet[jE]);
				if (sharedV != NULL && EdgeSet[iE].usedtimes == 1 && EdgeSet[jE].usedtimes == 1)
				{
					vector<Sensor*>assumedTriVertex;
					assumedTriVertex.push_back(EdgeSet[iE].start);
					assumedTriVertex.push_back(EdgeSet[iE].end);

					assumedTriVertex.push_back(EdgeSet[jE].start);
					assumedTriVertex.push_back(EdgeSet[jE].end);


					assumedTriVertex.erase(remove(assumedTriVertex.begin(),assumedTriVertex.end(),sharedV),assumedTriVertex.end());
					assert(assumedTriVertex.size() == 2);// size should be 2

					int linPos = findEdge(assumedTriVertex[0],assumedTriVertex[1]);
					if (linPos > -1 && EdgeSet[linPos].usedtimes == 1)
					{
						assert(EdgeSet[linPos].attachedTriangles.size() == 1);
						int triPos = findTriangle(sharedV,assumedTriVertex[0],assumedTriVertex[1]);
						if (triPos == -1)
						{
							++ n3Edges1usedtimeHole;
						}
						Triangle TTT(-1,sharedV,assumedTriVertex[0],assumedTriVertex[1]);
						auto it_removed = find_if(removedTris.begin(),removedTris.end(),[&](Triangle& tInS)
						{
							return (tInS == TTT);
						});
						if (it_removed!=removedTris.end())
						{
							cerr<<"Yes...";
						}
						
					}
				}
			}
		}
	}
						
	return n3Edges1usedtimeHole;
}
void Field::removeImpossibleHoles()
{
	for (int iE = 0; iE < EdgeSet.size(); ++ iE)
	{
		for ( int jE = iE + 1; jE < EdgeSet.size(); ++ jE )
		{
			if (EdgeSet[iE].TheSameLine(&(EdgeSet[jE])))
			{
				// the two line should not be the same line
				cerr<<"Same lines in the EdgeSet."<<endl;
				cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<endl;
				exit(EXIT_FAILURE);
			}
			else
			{
				//the two edges should sharing a vertex and the usedTimes have special relationship.
				Sensor* sharedV = EdgeSet[iE].sharedVertice(EdgeSet[jE]);
				if (sharedV != NULL && EdgeSet[iE].usedtimes == 1 && EdgeSet[jE].usedtimes == 1)
				{
					vector<Sensor*>assumedTriVertex;
					assumedTriVertex.push_back(EdgeSet[iE].start);
					assumedTriVertex.push_back(EdgeSet[iE].end);

					assumedTriVertex.push_back(EdgeSet[jE].start);
					assumedTriVertex.push_back(EdgeSet[jE].end);

					 
					assumedTriVertex.erase(remove(assumedTriVertex.begin(),assumedTriVertex.end(),sharedV),assumedTriVertex.end());
					assert(assumedTriVertex.size() == 2);// size should be 2
					
					int linPos = findEdge(assumedTriVertex[0],assumedTriVertex[1]);
					if (linPos > -1 && EdgeSet[linPos].usedtimes == 2)// find the wrong edge.
					{
						assert(EdgeSet[linPos].attachedTriangles.size() == 2);
						auto it_leastImportanceTri = min_element(EdgeSet[linPos].attachedTriangles.begin(),EdgeSet[linPos].attachedTriangles.end(),[&](Triangle* t1, Triangle* t2) -> bool
						{
							return (t1->triImportance < t2->triImportance);
						});

						if (it_leastImportanceTri!= EdgeSet[linPos].attachedTriangles.end())
						{
							Triangle* pTri = *it_leastImportanceTri;
							int edgePos12 = findEdge(pTri->t1,pTri->t2);
							int edgePos23 = findEdge(pTri->t2,pTri->t3);
							int edgePos13 = findEdge(pTri->t1,pTri->t3);

							assert(edgePos12 > -1 && edgePos23 > -1 && edgePos13 > -1);
							EdgeSet[edgePos12].usedtimes--;
							EdgeSet[edgePos23].usedtimes--;
							EdgeSet[edgePos13].usedtimes--;

							assert(EdgeSet[edgePos12].usedtimes >= 0);
							assert(EdgeSet[edgePos23].usedtimes >= 0);
							assert(EdgeSet[edgePos13].usedtimes >= 0);

							removedTris.push_back(*pTri);
							TriangleSet.erase(remove(TriangleSet.begin(),TriangleSet.end(),*pTri),TriangleSet.end());
							attachEdgesToEachLandmark();
							attachTriangleToEachLandmark();
							attachEdgesToTriangles();
							
						}
						else
						{
							cerr<<"Error at "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
							exit(EXIT_FAILURE);
						}
					}




				}


			}

		}
	}
}
void Field::removeItypeTeras(const int itype)
{
	for (int i = 0; i < nSensors; ++i)
	{
		Sensor* iSensor = &sensorPool[i];
		// Only landmarks have attached edges.
		if (iSensor->iType2Landmark[itype] && iSensor->iTypeAttachedTriangles[itype].size() == 3)
		{
			// deal with the end vertices's.
			set<Sensor*>endVertex;
			for (auto it_iTri = iSensor->iTypeAttachedTriangles[itype].begin();
				it_iTri != iSensor->iTypeAttachedTriangles[itype].end(); ++ it_iTri)
			{
				Triangle* iTri = *it_iTri;
				if (iTri->t1 != iSensor)
				{
					endVertex.insert(iTri->t1);
				}
				if (iTri->t2 != iSensor)
				{
					endVertex.insert(iTri->t2);
				}
				if (iTri->t3 != iSensor)
				{
					endVertex.insert(iTri->t3);
				}

			}

			// if the size of endVetex is 3, then there is a tera.
			// remove one face of the tera with smallest importance.
			if (endVertex.size() == 3)
			{
				auto it_end_fourth_tri = endVertex.begin();
				Sensor* fourthTri_v1 = *(it_end_fourth_tri);
				Sensor* fourthTri_v2 = *(++it_end_fourth_tri);
				Sensor* fourthTri_v3 = *(++it_end_fourth_tri);

				vector<Triangle> teraTriFaces;
				// Tri_face1234
				Triangle iS_v1_v2(-1,iSensor,fourthTri_v1,fourthTri_v2);
				Triangle iS_v2_v3(-1,iSensor,fourthTri_v2,fourthTri_v3);
				Triangle iS_v1_v3(-1,iSensor,fourthTri_v1,fourthTri_v3);
				Triangle v1_v2_v3(-1,fourthTri_v1,fourthTri_v2,fourthTri_v3);

				teraTriFaces.push_back(iS_v1_v2);
				teraTriFaces.push_back(iS_v2_v3);
				teraTriFaces.push_back(iS_v1_v3);
				teraTriFaces.push_back(v1_v2_v3);

				auto it_leastImportanceTri = min_element(teraTriFaces.begin(),teraTriFaces.end(),[&](Triangle& t1, Triangle& t2) -> bool
				{
					return (t1.triImportance < t2.triImportance);
				});
				// Find the least importance triangle in the triangleSet.
				auto it_inTriSet = find_if(iTypeTriangleSet[itype].begin(),iTypeTriangleSet[itype].end(),[&](Triangle& tri)->bool
				{
					return (tri == *it_leastImportanceTri);
				});
				// must find it, or it is a wrong triangle.
				if (it_inTriSet != iTypeTriangleSet[itype].end())
				{
					// Update the used times of the triangle edges.
					// tmp line segment to -- usedtimes.
					LineSegment L12(-1,(*it_inTriSet).t1,(*it_inTriSet).t2);
					LineSegment L23(-1,(*it_inTriSet).t2,(*it_inTriSet).t3);
					LineSegment L13(-1,(*it_inTriSet).t1,(*it_inTriSet).t3);

					auto it_edge12 = find_if(iTypeEdgesSet[itype].begin(),iTypeEdgesSet[itype].end(),[&](LineSegment& L)
					{
						return (L == L12); 		
					});
					auto it_edge23 = find_if(iTypeEdgesSet[itype].begin(),iTypeEdgesSet[itype].end(),[&](LineSegment& L)
					{
						return (L == L23); 		
					});
					auto it_edge13 = find_if(iTypeEdgesSet[itype].begin(),iTypeEdgesSet[itype].end(),[&](LineSegment& L)
					{
						return (L == L13); 		
					});

					assert((*it_edge12).usedtimes == 2);
					assert((*it_edge23).usedtimes == 2);
					assert((*it_edge13).usedtimes == 2);

					(*it_edge12).usedtimes--;
					(*it_edge23).usedtimes--;
					(*it_edge13).usedtimes--;

					// Remove the triangle
					iTypeRemovedTris[itype].push_back(*it_inTriSet);
					TriangleSet.erase(remove(iTypeTriangleSet[itype].begin(),iTypeTriangleSet[itype].end(),*it_inTriSet),iTypeTriangleSet[itype].end());
					
					itypeAttachEdgesToEachLandmark(itype);
					itypeAttachTriangleToEachLandmark(itype);
					itypeAttachEdgesToTriangles(itype);

				}
				
			}

		}
	}
}
void Field::removeTeras()
{
	// Removing tera by delete a triangle.
	for (int i = 0; i < nSensors; ++i)
	{
		Sensor* iSensor = &sensorPool[i];
		// Only landmarks have attached edges.
		if (iSensor->landmark && iSensor->attachedTriangles.size() == 3)
		{
			// deal with the end vertices's.
			set<Sensor*>endVertex;
			for (auto it_iTri = iSensor->attachedTriangles.begin();
				it_iTri != iSensor->attachedTriangles.end(); ++ it_iTri)
			{
				Triangle* iTri = *it_iTri;
				if (iTri->t1 != iSensor)
				{
					endVertex.insert(iTri->t1);
				}
				if (iTri->t2 != iSensor)
				{
					endVertex.insert(iTri->t2);
				}
				if (iTri->t3 != iSensor)
				{
					endVertex.insert(iTri->t3);
				}

			}

			// if the size of endVetex is 3, then there is a tera.
			// remove one face of the tera with smallest importance.
			if (endVertex.size() == 3)
			{
				auto it_end_fourth_tri = endVertex.begin();
				Sensor* fourthTri_v1 = *(it_end_fourth_tri);
				Sensor* fourthTri_v2 = *(++it_end_fourth_tri);
				Sensor* fourthTri_v3 = *(++it_end_fourth_tri);

				vector<Triangle> teraTriFaces;
				// Tri_face1234
				Triangle iS_v1_v2(-1,iSensor,fourthTri_v1,fourthTri_v2);
				Triangle iS_v2_v3(-1,iSensor,fourthTri_v2,fourthTri_v3);
				Triangle iS_v1_v3(-1,iSensor,fourthTri_v1,fourthTri_v3);
				Triangle v1_v2_v3(-1,fourthTri_v1,fourthTri_v2,fourthTri_v3);

				teraTriFaces.push_back(iS_v1_v2);
				teraTriFaces.push_back(iS_v2_v3);
				teraTriFaces.push_back(iS_v1_v3);
				teraTriFaces.push_back(v1_v2_v3);

				auto it_leastImportanceTri = min_element(teraTriFaces.begin(),teraTriFaces.end(),[&](Triangle& t1, Triangle& t2) -> bool
				{
					return (t1.triImportance < t2.triImportance);
				});
				// Find the least importance triangle in the triangleSet.
				auto it_inTriSet = find_if(TriangleSet.begin(),TriangleSet.end(),[&](Triangle& tri)->bool
				{
					return (tri == *it_leastImportanceTri);
				});
				// must find it, or it is a wrong triangle.
				if (it_inTriSet != TriangleSet.end())
				{
					// Update the used times of the triangle edges.
					// tmp linesegment to -- usedtimes.
					LineSegment L12(-1,(*it_inTriSet).t1,(*it_inTriSet).t2);
					LineSegment L23(-1,(*it_inTriSet).t2,(*it_inTriSet).t3);
					LineSegment L13(-1,(*it_inTriSet).t1,(*it_inTriSet).t3);

					auto it_edge12 = find_if(EdgeSet.begin(),EdgeSet.end(),[&](LineSegment& L)
					{
						return (L == L12); 		
					});
					auto it_edge23 = find_if(EdgeSet.begin(),EdgeSet.end(),[&](LineSegment& L)
					{
						return (L == L23); 		
					});
					auto it_edge13 = find_if(EdgeSet.begin(),EdgeSet.end(),[&](LineSegment& L)
					{
						return (L == L13); 		
					});

					assert((*it_edge12).usedtimes == 2);
					assert((*it_edge23).usedtimes == 2);
					assert((*it_edge13).usedtimes == 2);

					(*it_edge12).usedtimes--;
					(*it_edge23).usedtimes--;
					(*it_edge13).usedtimes--;

					// Remove the triangle
					removedTris.push_back(*it_inTriSet);
					TriangleSet.erase(remove(TriangleSet.begin(),TriangleSet.end(),*it_inTriSet),TriangleSet.end());
					attachEdgesToEachLandmark();
					attachTriangleToEachLandmark();
					attachEdgesToTriangles();
				}
//  				else
//  				{
//  					cerr<<"Tera exist, but can not delete a face"<<endl;
//  					exit(EXIT_FAILURE);
//  				}
			}

		}
	}
}
void Field::GenerateVoronoiCellForSkeleton()
{
	int i, j;
	int count = 0;
	for(i = 0; i < nSensors; ++i)
	{
		Sensor *b = &sensorPool[i];
		if(!b->landmark && b->onEdge_Degree == 0)
		{
			if (b->index == 58208)
			{
				printf("wait\n");
			}
			b->LocalFloodingForEstablishCellforSkeleton(1);
			if(b->VoronoiCellID.size() > 1)
			{
				++count;
				VoronoiCellSensorSet.insert(b);
				printf("VoronoiCellSensorSet %d\n", b->index);
			}
		}
	}
	printf("Voronoi Cell Set for skeleton %d\n", count);
	printf("Voronoi Cell for skeleton has defined!\n");
	//Get landmark neighbors for each landmark
	SCITER VoronoiCellEdgeSetpos;
	set<int>::iterator VoronoiCellIDpos;
	for(VoronoiCellEdgeSetpos = VoronoiCellSensorSet.begin(); VoronoiCellEdgeSetpos != VoronoiCellSensorSet.end(); ++VoronoiCellEdgeSetpos)
	{
		Sensor *s = *VoronoiCellEdgeSetpos;
		int SIZE = s->VoronoiCellID.size();

		int ID[128] = {0};
		for(i = 0, VoronoiCellIDpos = s->VoronoiCellID.begin(); VoronoiCellIDpos != s->VoronoiCellID.end(); ++VoronoiCellIDpos, ++i)
		{
			ID[i] = *VoronoiCellIDpos;
		}
		for(i = 0; i < SIZE; ++i)
		{
			for(j = 0; j < SIZE; ++j)
			{
				if(ID[i] != ID[j])
				{
					sensorPool[ID[i]].LandmarkNeighbors.insert(&sensorPool[ID[j]]);
				}
			}
		}
	}
}

#if 0
void Field::RefineSurface()
{
	GetRefineSensorAndEdgeSet();
	bool updatetag = false;
	SCITER RefineSensorSetpos2;
	SCITER RefineSensorSetpos1;
	int LineIDindex = EdgeSet.size();
	set<LineSegment*> EnclosedLines;
	set<LineSegment*>::iterator EnclosedLinespos;
	do 
	{
		updatetag = false;
		EnclosedLines.clear();
		for(RefineSensorSetpos1 = RefineSensorSet.begin(); RefineSensorSetpos1 != RefineSensorSet.end(); ++RefineSensorSetpos1)
		{
			Sensor *R = *RefineSensorSetpos1;
			for(RefineSensorSetpos2 = RefineSensorSet.begin(); RefineSensorSetpos2 != RefineSensorSet.end(); ++RefineSensorSetpos2)
			{
				bool Lcanceled = false;
				Sensor *l = *RefineSensorSetpos2;
				if(/*LandmarkBelongToRefineSensorSet(l) && */R->index != l->index)
				{
					LineSegment NewL(LineIDindex, R, l);
					for(unsigned int j = 0; j < EdgeSet.size(); ++j)
					{
						if(NewL.TheSameLine(&EdgeSet[j]))
						{
							Lcanceled= true;
							break;
						}	
					}
					if(!Lcanceled)
					{
						EnclosedLines =  TriangleEnclosed(&NewL);
						if(!EnclosedLines.empty())
						{
							for(EnclosedLinespos = EnclosedLines.begin(); EnclosedLinespos != EnclosedLines.end(); ++EnclosedLinespos)
							{
								LineSegment *l = *EnclosedLinespos;
								l->usedtimes = l->The3rdPoint.size();
							}
							EdgeSet.push_back(NewL);
							if(NewL.usedtimes < 2)
							{
								RefineEdgeSet.push_back(&NewL);
							}
							else
								printf("NewL used 4 times directly\n");
							UpdateRefineEdgeSet();
							updatetag = true;
							break;
						}
					}
				}
			}
			if(updatetag)
				break;
		}
		printf("RefineEdgeSetSize is %d\n", RefineEdgeSet.size());
	} while (!RefineEdgeSet.empty());
}
#endif

void Field::RefineSurface()
{
	GetRefineEdgeSet();
	if(!RefineEdgeSet.empty())
	{
		DetectHoles();
		ClassifyHoleSensors();
		FillHoles();
		//GetRefineEdgeSet();
	}
	
	//updating triangle neighborhood
	for (unsigned int i = 0; i < TriangleSet.size(); ++i)
	{
		TriangleSet[i].EstablishNeighborhood();
	}
}

void Field::RefineSurfaceForSkeleton()
{
	int debugnum = 0;
	vector<Sensor*>::iterator LandmarkSetpos;
	bool CancelTheLine = false;
	for(LandmarkSetpos = LandmarkSet.begin(); LandmarkSetpos != LandmarkSet.end(); ++LandmarkSetpos)
	{
		Sensor *CurrentLandmark = *LandmarkSetpos;
		
		for(SCITER LNpos = CurrentLandmark->LandmarkNeighbors.begin(); LNpos != CurrentLandmark->LandmarkNeighbors.end(); ++LNpos)
		{
			CancelTheLine = false;
			Sensor *LN = *LNpos;
			if (CurrentLandmark->NotConnectWith(LN))
			{
				Path line(CurrentLandmark, LN);
				line.EstablishDirectPath(CurrentLandmark, LN);
				for(unsigned int p = 1; p < line.path.size()-1; ++p)
				{
					Sensor *lp = line.path[p];
					if(lp->BelongOtherIntersectLine(line.start->index, line.end->index))
					{
						CancelTheLine = true;
						printf("Detect crossing!\n");
						break;
					}


					for(SCITER lpnpos = lp->neighbors.begin(); lpnpos != lp->neighbors.end(); ++lpnpos)
					{
						Sensor *lpn = *lpnpos;
						if(lpn->onEdge_Degree != 0 && lpn->BelongOtherIntersectLine(line.start->index, line.end->index))
						{
							CancelTheLine = true;
							printf("Detect crossing!\n");
							break;
						}
					}

				}
				if(!CancelTheLine)
				{
					++debugnum;
					if(debugnum == 17)
						printf("wait\n");
					CurrentLandmark->ConnectedID.insert(LN->index);
					LN->ConnectedID.insert(CurrentLandmark->index);
					FILE *fp;
					errno_t err;
					err = fopen_s(&fp, "Line", "a");
					if (!err) 
					{
						fprintf_s(fp, "$");
					}
					else 
					{
						fprintf_s(stderr, "open file %s to write failed!\n", "Line");
						return;
					}
					fclose(fp);
					dumpLine("Line", CurrentLandmark, LN);
					for(unsigned int p = 1; p < line.path.size()-1; ++p)
					{
						Sensor *lp = line.path[p];
						lp->line_ends_ID.insert(make_pair(CurrentLandmark->index, LN->index));
					}
				}
			}
		}
	}
}

void Field::FillHoles()
{
	SCITER HoleSensorSetPos;
	Sensor* cutpair[2] = {NULL, NULL};
	//Sensor* original_cutpair[2] = {NULL, NULL};
	//vector<Triangle> TempTriangleSet;
	unsigned int cutpair_index;
	unsigned int triangle_sensor_count;
	unsigned int i, j;
	int Tri_ID = TriangleSet.size();
	bool FillWayChange = false;

	for ( i = 0; i < HoleSet.size(); ++i )
	{
		/*
		for(HoleSensorSetPos = HoleSet[i].HoleSensorSet.begin(); HoleSensorSetPos != HoleSet[i].HoleSensorSet.end(); ++HoleSensorSetPos)
		{
		Sensor *temps = *HoleSensorSetPos;
		temps->updateHoleSensorNeighbors(HoleSet[i].HoleSensorSet);
		}*/
		//TempTriangleSet.clear();
		triangle_sensor_count = 2;
		cutpair_index = 0;
		if(HoleSet[i].HoleSensorSet.size() > 10)
		{
			UniformRandVar Urand(0, HoleSet[i].HoleSensorSet.size()-1);
			int iter_count = (int)Urand.value();
			HoleSensorSetPos = HoleSet[i].HoleSensorSet.begin();
			while (iter_count)
			{
				++HoleSensorSetPos;
				--iter_count;
			}
		}
		else
			HoleSensorSetPos = HoleSet[i].HoleSensorSet.begin();
		Sensor *s = *HoleSensorSetPos;
		if(i == 32)
			bool stop  = 1;
		s->HoleSensorFlooding(HoleSet[i].HoleSensorSet);
		//////////////////////////////////////////////////////////////////////////
		//get cut pair
		for ( HoleSensorSetPos = HoleSet[i].HoleSensorSet.begin(); HoleSensorSetPos != HoleSet[i].HoleSensorSet.end(); ++HoleSensorSetPos )
		{
			Sensor *sen = *HoleSensorSetPos;
			if ( sen->children.size() == 0 )
			{
				cutpair[cutpair_index] = sen;
				++cutpair_index;
			}
		}
		
		//////////////////////////////////////////////////////////////////////////
		while ( triangle_sensor_count != HoleSet[i].HoleSensorSet.size() )
		{
			FillWayChange = false;
			int alternative_parent = 0xffff;
			int current_parent = 0xffff;
			int final_parent_cut_index = 0xffff;
			//LineSegment l1(0xffff, cutpair[0]->parent, cutpair[0]);
			//LineSegment l2(0xffff, cutpair[0]->parent, cutpair[1]);
			LineSegment l1(0xffff, cutpair[0], cutpair[1]);
			LineSegment l2;
			LineSegment l3;
			/*
			LineSegment *targetedge1 = NULL;
			LineSegment *targetedge2 = NULL;
			LineSegment *targetedge3 = NULL;
			*/
			Path p1;
			Path p2;
			if(cutpair[0]->hole_parent != NULL)
				p1.EstablishDirectPath(cutpair[0]->hole_parent, cutpair[1]);
			if(cutpair[1]->hole_parent != NULL)
				p2.EstablishDirectPath(cutpair[1]->hole_parent, cutpair[0]);
			if(!p1.path.empty() && !p2.path.empty())
			{
				if(p1.path.size() > p2.path.size())
				{
					l2.init(0xffff, cutpair[1]->hole_parent, cutpair[0]);
					l3.init(0xffff, cutpair[1]->hole_parent, cutpair[1]);
					current_parent = 1;
					alternative_parent = 0;
				}
				else
				{
					l2.init(0xffff, cutpair[0]->hole_parent, cutpair[0]);
					l3.init(0xffff, cutpair[0]->hole_parent, cutpair[1]);
					current_parent = 0;
					alternative_parent = 1;
				}
			}
			else
			{
				if(!p1.path.empty())
				{
					l2.init(0xffff, cutpair[0]->hole_parent, cutpair[0]);
					l3.init(0xffff, cutpair[0]->hole_parent, cutpair[1]);
					current_parent = 0;
					alternative_parent = 1;
				}
				else
				{
					l2.init(0xffff, cutpair[1]->hole_parent, cutpair[0]);
					l3.init(0xffff, cutpair[1]->hole_parent, cutpair[1]);
					current_parent = 1;
					alternative_parent = 0;
				}
			}
			for(j = 0; j < EdgeSet.size(); ++j)
			{
				if( (l1.TheSameLine(&EdgeSet[j])) /*|| l1.Lpath.PathCrossing(1)*/)
				{
					//targetedge1 = &EdgeSet[j];
					if(EdgeSet[j].usedtimes == 2)
					{
						FillWayChange = true;
						break;
					}
				}
				if( (l2.TheSameLine(&EdgeSet[j])) /*|| l2.Lpath.PathCrossing(1)*/)
				{
					//targetedge2 = &EdgeSet[j];
					if(EdgeSet[j].usedtimes == 2)
					{
						FillWayChange = true;
						break;
					}
				}
				if( (l3.TheSameLine(&EdgeSet[j])) /*|| l3.Lpath.PathCrossing(1)*/)
				{
					//targetedge3 = &EdgeSet[j];
					if(EdgeSet[j].usedtimes == 2)
					{
						FillWayChange = true;
						break;
					}
				}
			}
			if(!FillWayChange || cutpair[alternative_parent]->hole_parent == NULL)
			{
				Triangle T1(Tri_ID, cutpair[current_parent]->hole_parent, cutpair[0], cutpair[1]);
				/*
				if(targetedge1 != NULL)
				{
					targetedge1->The3rdPoint.push_back(cutpair[current_parent]->hole_parent);
					targetedge1->usedtimes = targetedge1->The3rdPoint.size();
				}
				else
				{
					l1.The3rdPoint.push_back(cutpair[current_parent]->hole_parent);
					l1.usedtimes = l1.The3rdPoint.size();
					EdgeSet.push_back(l1);
				}
				if(targetedge2 != NULL)
				{
					targetedge2->The3rdPoint.push_back(cutpair[1]);
					targetedge2->usedtimes = targetedge2->The3rdPoint.size();
				}
				else
				{
					l2.The3rdPoint.push_back(cutpair[1]);
					l2.usedtimes = l2.The3rdPoint.size();
					EdgeSet.push_back(l2);
				}
				if(targetedge3 != NULL)
				{
					targetedge3->The3rdPoint.push_back(cutpair[0]);
					targetedge3->usedtimes = targetedge3->The3rdPoint.size();
				}
				else
				{
					l3.The3rdPoint.push_back(cutpair[0]);
					l3.usedtimes = l3.The3rdPoint.size();
					EdgeSet.push_back(l3);
				}
				*/
				final_parent_cut_index = current_parent;
				T1.Tri_field = this;
				T1.GetNormalVector();
				TriangleSet.push_back(T1);
				//dumpTriangle("triangle", &T1);
				++Tri_ID;
				++Total_Number_of_Triangles;
				printf("Generate %d triangles\n", Total_Number_of_Triangles);
				if(Total_Number_of_Triangles == 115)
					bool stop = 1;
				++triangle_sensor_count;
			}
			else
			{
				if(cutpair[alternative_parent]->hole_parent == NULL)
				{
					++triangle_sensor_count;
					break;
				}
				Triangle T1(Tri_ID, cutpair[alternative_parent]->hole_parent, cutpair[0], cutpair[1]);
				/*
				if(targetedge1 != NULL)
				{
					targetedge1->The3rdPoint.push_back(cutpair[alternative_parent]->hole_parent);
					targetedge1->usedtimes = targetedge1->The3rdPoint.size();
				}
				else
				{
					l1.The3rdPoint.push_back(cutpair[alternative_parent]->hole_parent);
					l1.usedtimes = l1.The3rdPoint.size();
					EdgeSet.push_back(l1);
				}
				if(targetedge2 != NULL)
				{
					targetedge2->The3rdPoint.push_back(cutpair[1]);
					targetedge2->usedtimes = targetedge2->The3rdPoint.size();
				}
				else
				{
					l2.The3rdPoint.push_back(cutpair[1]);
					l2.usedtimes = l2.The3rdPoint.size();
					EdgeSet.push_back(l2);
				}
				if(targetedge3 != NULL)
				{
					targetedge3->The3rdPoint.push_back(cutpair[0]);
					targetedge3->usedtimes = targetedge3->The3rdPoint.size();
				}
				else
				{
					l3.The3rdPoint.push_back(cutpair[0]);
					l3.usedtimes = l3.The3rdPoint.size();
					EdgeSet.push_back(l3);
				}
				*/
				final_parent_cut_index = alternative_parent;
				T1.Tri_field = this;
				T1.GetNormalVector();
				TriangleSet.push_back(T1);
				//dumpTriangle("triangle", &T1);
				++Tri_ID;
				++Total_Number_of_Triangles;
				if(Total_Number_of_Triangles == 30)
					printf("wait\n");
				printf("Generate %d triangles\n", Total_Number_of_Triangles);
				if(Total_Number_of_Triangles == 115)
					bool stop = 1;
				++triangle_sensor_count;
			}
			Sensor *cut0 = cutpair[final_parent_cut_index]->hole_parent;
			Sensor *cut1 = cutpair[1-final_parent_cut_index];
			cutpair[0] = cut0;
			cutpair[1] = cut1;
		}
		for ( HoleSensorSetPos = HoleSet[i].HoleSensorSet.begin(); HoleSensorSetPos != HoleSet[i].HoleSensorSet.end(); ++HoleSensorSetPos )
		{
			Sensor *sen = *HoleSensorSetPos;
			sen->candidateParents.clear();
		}
	}
}

void Field::GetRefineEdgeSet()
{
	HoleSet.clear();
	for(unsigned int i = 0; i < EdgeSet.size(); ++i)
	{
		if(EdgeSet[i].usedtimes < 2)
		{
			RefineEdgeSet.push_back(&EdgeSet[i]);
		}
		if(EdgeSet[i].usedtimes < 1)
		{
			printf("There are problems in edge %d\n", i);
			exit(0);
		}
	}
}

void Field::DetectHoles()
{
	set<Sensor*> hole_sensor_set;
	vector<LineSegment*> hole_edge_set;				//used to store confirmed hole edge
	vector<LineSegment*> candidate_hole_edge_set;	//used to store probable hole edge
	vector<LineSegment*> LineBranch;
	vector<int> LineBranchIndex;
	Sensor *Joint = NULL;
	Sensor *PreviousJoint = NULL;

	int next_start_times = 0;
	int next_end_times = 0;
	int HoleCount = 0;
	int forbidden_index = 0xffff;
	int end_index = 0xffff;

	LineSegment *lstart = NULL;
	LineSegment *lcurrent = NULL;
	LineSegment *lnext = NULL;
	unsigned int i;
	unsigned int j;

	do 
	{
		bool continue_tag = false;
		lstart = RefineEdgeSet[0];
		forbidden_index = 0xffff;
		end_index = lstart->end->index;
		lcurrent = lstart;
		candidate_hole_edge_set.clear();
		LineBranchIndex.clear();
		hole_sensor_set.clear();
		hole_edge_set.clear();
		hole_edge_set.push_back(lcurrent);
		candidate_hole_edge_set.push_back(lstart);
		do
		{	
			continue_tag = false;
			LineBranch.clear();
			next_start_times = 0;
			next_end_times = 0;
			for(j = 0; j < RefineEdgeSet.size(); ++j)
			{
				if( (*lcurrent) == (*RefineEdgeSet[j]) )
					continue;
				if(!lcurrent->start->LineSegmentNeighborState && lcurrent->start->index != forbidden_index)
				{
					if(lcurrent->start->index == RefineEdgeSet[j]->start->index || lcurrent->start->index == RefineEdgeSet[j]->end->index)
					{
						if(LineBranch.empty())
						{
							PreviousJoint = Joint;
						}
						Joint = lcurrent->start;
						LineBranch.push_back(RefineEdgeSet[j]);
						++next_start_times;
					}

				}
				else if(!lcurrent->end->LineSegmentNeighborState && lcurrent->end->index != forbidden_index)
				{
					if(lcurrent->end->index == RefineEdgeSet[j]->start->index || lcurrent->end->index == RefineEdgeSet[j]->end->index)
					{
						if(LineBranch.empty())
						{
							PreviousJoint = Joint;
						}
						Joint = lcurrent->end;
						LineBranch.push_back(RefineEdgeSet[j]);
						++next_end_times;

					}

				}
				else
					continue;
			}
			if(PreviousJoint != NULL)
				PreviousJoint->LineSegmentNeighborState = false;
			if(LineBranch.size() > 1)
			{
				if(next_start_times > 1)
				{
					forbidden_index = lcurrent->start->index;
				}
				if(next_end_times > 1)
				{
					forbidden_index = lcurrent->end->index;
				}
				if(next_end_times > 1 && next_start_times > 1)
				{
					printf("next times wrong!\n");
					exit(0);
				}
				LineBranchIndex.push_back(forbidden_index);
				bool END = false;
				for(i = 0; i < LineBranch.size(); ++i)
				{
					if(LineBranch[i]->start->index == end_index || LineBranch[i]->end->index == end_index)
					{
						forbidden_index = 0xffff;
						if(Joint != NULL)
						{
							Joint->LineSegmentNeighborState = true;
						}
						lnext = LineBranch[i];
						lcurrent = lnext;
						hole_edge_set.push_back(lcurrent);
						candidate_hole_edge_set.push_back(lcurrent);
						END = true;
						break;
					}
				}
				if(!END)
				{
					lstart = LineBranch[0];
					lstart->start->LineSegmentNeighborState = false;
					lstart->end->LineSegmentNeighborState = false;
					end_index = forbidden_index;
					lcurrent = lstart;
					hole_edge_set.clear();
					hole_edge_set.push_back(lcurrent);
					candidate_hole_edge_set.push_back(lstart);
					continue_tag = true;
				}
				//bool breaktag = false;
				int LBIcount = 0;
				for(i = 0; i < LineBranchIndex.size(); ++i)
				{
					if(forbidden_index == LineBranchIndex[i])
					{
						++LBIcount;
					}
				}
				if(LBIcount == 2)
				{
					int start_index = 0;
					hole_edge_set.clear();
					for(j = 0; j < candidate_hole_edge_set.size(); ++j)
					{
						if(candidate_hole_edge_set[j]->start->index == forbidden_index || candidate_hole_edge_set[j]->end->index == forbidden_index)
						{
							start_index = j+1;
							break;
						}
					}
					for(j = start_index; j < candidate_hole_edge_set.size()-1; ++j)
					{
						hole_edge_set.push_back(candidate_hole_edge_set[j]);
					}

					break;
				}
			}
			//choose lnext
			else
			{
				forbidden_index = 0xffff;
				if(Joint != NULL)
				{
					Joint->LineSegmentNeighborState = true;
				}
				lnext = LineBranch[0];
				lcurrent = lnext;
				hole_edge_set.push_back(lcurrent);
				candidate_hole_edge_set.push_back(lcurrent);
			}
		} while( (end_index != lcurrent->start->index && end_index != lcurrent->end->index) || continue_tag);

		//Get rid of the last edge in hole_edge_set, because it is the same as the start line
		//vector<LineSegment*>::iterator hole_edge_set_pos;
		//hole_edge_set_pos = hole_edge_set.end();
		//--hole_edge_set_pos;
		//hole_edge_set.erase(hole_edge_set_pos);
		//////////////////////////////////////////////////////////////////////////

		for(i = 0; i < hole_edge_set.size(); ++i)
		{
			hole_sensor_set.insert(hole_edge_set[i]->start);
			if(hole_edge_set[i]->start->LineSegmentNeighborState)
			{
				hole_edge_set[i]->start->LineSegmentNeighborState = false;
			}
			hole_sensor_set.insert(hole_edge_set[i]->end);
			if(hole_edge_set[i]->end->LineSegmentNeighborState)
			{
				hole_edge_set[i]->end->LineSegmentNeighborState = false;
			}
			//Get HoleSetNeighbors

			hole_edge_set[i]->start->HoleSetNeighbors.insert(hole_edge_set[i]->end);
			hole_edge_set[i]->end->HoleSetNeighbors.insert(hole_edge_set[i]->start);
		}

		HoleSet.push_back(Hole(HoleCount, hole_sensor_set));

		printf("Detected a hole, HoleSet size is %d\n", HoleSet.size());
		++HoleCount;
		UpdateRefineEdgeSet(hole_edge_set);
		printf("The size of refining edge set is %d\n", RefineEdgeSet.size());
	} while (!RefineEdgeSet.empty());
}


void Field::GetRefineSensorAndEdgeSet()
{
	for(unsigned int i = 0; i < EdgeSet.size(); ++i)
	{
		if(EdgeSet[i].usedtimes < 2)
		{
			RefineEdgeSet.push_back(&EdgeSet[i]);
			RefineSensorSet.insert(EdgeSet[i].start);
			RefineSensorSet.insert(EdgeSet[i].end);
		}
		if(EdgeSet[i].usedtimes < 1)
		{
			printf("There are problems in edge %d\n", i);
			exit(0);
		}
	}
}

void Field::ClassifyHoleSensors()
{
	unsigned int i = 0;
	for(i = 0; i < HoleSet.size(); ++i)
	{
		SCITER HoleSensorPos;
		for(HoleSensorPos = HoleSet[i].HoleSensorSet.begin(); HoleSensorPos != HoleSet[i].HoleSensorSet.end(); ++HoleSensorPos)
		{
			Sensor *st = *HoleSensorPos;
			st->BelongInHoleID.push_back(HoleSet[i].holeID);
		}
	}
}

#if 0
void Field::UpdateRefineEdgeSet()
{
	vector<LineSegment*>::iterator RefineEdgeSetPos;

	do
	{
		for(RefineEdgeSetPos = RefineEdgeSet.begin(); RefineEdgeSetPos != RefineEdgeSet.end(); ++RefineEdgeSetPos)
		{
			LineSegment *l = *RefineEdgeSetPos;
			if(l->usedtimes == 2)
			{
				RefineEdgeSet.erase(RefineEdgeSetPos);
				break;
			}
		}
	} while(!EndUpdating());
}
#endif

void Field::UpdateRefineEdgeSet(vector<LineSegment*> HoleEdgeSet)
{
	unsigned int i;
	int endcount = 0;
	bool erasetag = false;
	vector<LineSegment*>::iterator RefineEdgeSetPos;
	do
	{
		erasetag = false;
		for(RefineEdgeSetPos = RefineEdgeSet.begin(); RefineEdgeSetPos != RefineEdgeSet.end(); ++RefineEdgeSetPos)
		{
			LineSegment *l = *RefineEdgeSetPos;
			for(i = 0; i < HoleEdgeSet.size(); ++i)
			{
				if(l->LineID == HoleEdgeSet[i]->LineID)
				{
					RefineEdgeSet.erase(RefineEdgeSetPos);
					++endcount;
					erasetag = true;
					break;
				}
			}
			if(erasetag)
				break;
		}
	} while(endcount != HoleEdgeSet.size());

}

set<LineSegment*> Field::TriangleEnclosed(LineSegment *l)
{
	set<LineSegment*> enclosed;
	set<LineSegment*>::iterator enclosedpos;
	vector<LineSegment*> lstartEdgeSet;
	vector<LineSegment*> lendEdgeSet;
	//search enclosed lines
	unsigned int i, j;
	for(i = 0; i < RefineEdgeSet.size(); ++i)
	{
		//looking for l->start point's neighbor edge
		if(RefineEdgeSet[i]->start->index == l->start->index || RefineEdgeSet[i]->end->index == l->start->index)
		{
			lstartEdgeSet.push_back(RefineEdgeSet[i]);
		}
		//looking for l->end point's neighbor edge
		if(RefineEdgeSet[i]->start->index == l->end->index || RefineEdgeSet[i]->end->index == l->end->index)
		{
			lendEdgeSet.push_back(RefineEdgeSet[i]);
		}
	}
 	for(i = 0; i < lstartEdgeSet.size(); ++i)
	{
		for(j = 0; j < lendEdgeSet.size(); ++j)
		{
			if(lstartEdgeSet[i]->start->index == lendEdgeSet[j]->start->index 
				|| lstartEdgeSet[i]->start->index == lendEdgeSet[j]->end->index 
				|| lstartEdgeSet[i]->end->index == lendEdgeSet[j]->start->index 
				|| lstartEdgeSet[i]->end->index == lendEdgeSet[j]->end->index)
			{
				enclosed.insert(lstartEdgeSet[i]);
				enclosed.insert(lendEdgeSet[j]);
				if(i == 0 && enclosed.size() == 2)
					break;
			}
		}
		if(i == 0 && enclosed.size() == 2)
			break;
	}
	
	//Add the 3rd point of enclosed lines & l
	if(!enclosed.empty())
	{
		//Add the 3rd point of l
		enclosedpos = enclosed.begin();
		LineSegment *temp = *enclosedpos;
		if(temp->start->index != l->start->index && temp->start->index != l->end->index)
			l->The3rdPoint.push_back(temp->start);
		else
			l->The3rdPoint.push_back(temp->end);
		l->usedtimes = l->The3rdPoint.size();
		for(enclosedpos = enclosed.begin(); enclosedpos != enclosed.end(); ++enclosedpos)
		{
			LineSegment *en = *enclosedpos;
			if(l->start->index != en->start->index && l->start->index != en->end->index)
			{
				en->The3rdPoint.push_back(l->start);
			}
			if(l->end->index != en->start->index && l->end->index != en->end->index)
			{
				en->The3rdPoint.push_back(l->end);
			}
		}
		
	}
	return enclosed;
}

void Field::dumpTopology(char * filename)			
{
	if (!filename)
		return;

	FILE *fp;
	errno_t err;
	err = fopen_s(&fp, filename, "w");
	if (!err)
	{
		//nsensors
		fprintf_s(fp, "#%d %lf %lf\n", globalField->nSensors, avg1HopNeighborhoodSize, avgNHopNeighborhoodSize);
		for(int i = 0; i < globalField->nSensors; i++) 
		{
			Sensor &s = globalField->sensorPool[i];
			fprintf_s(fp, "%d %lf %lf %d %d %d %d %d %d %d %lf %lf %lf %d\n", i, s.onEdge_Degree, s.Criticality, s.landmark, s.convex, s.concave, s.saddle, s.asSkeleton, s.Out, s.VoronoiCellID.size(), s.location.x, s.location.y, s.location.z,s.boundarySensorType);
		}
	}
	else 
	{
		fprintf_s(stderr, "open file %s to write failed!\n", filename);
		return;
	}
	fclose(fp);
}

void Field::dumpTopologyFull(char * filename)
{

	if (!filename) 
		return;

	FILE *fp;
	errno_t err;
	err = fopen_s(&fp, filename, "w");

	if (!err) 
	{
		for(int i = 0; i < globalField->nSensors; i++) 
		{
			Sensor &s = globalField->sensorPool[i];

			fprintf_s(fp, "%d %lf %lf %d %d %d %d\n", i, s.onEdge_Degree, 
				s.Criticality, s.concave ? 1 : 0,
				s.convex ? 1 : 0, s.saddle ? 1 : 0, s.asSkeleton ? 1:0);

		}
	}
	else 
	{
		fprintf_s(stderr, "open file %s to write failed!\n", filename);
		return;
	}
	fclose(fp);
}

void Field::recoverPreprocessedNetwork(char* topoFile)
{
	if(topoFile) {
		FILE *fp;
		errno_t err;
		err = fopen_s(&fp, topoFile, "r");
		char line[1024] = {0};
		if (err)  {
			fprintf_s(stderr, "no topo file %s found\n", topoFile);
			exit(0);
		}
		do {
			memset(line, 0, sizeof(line));
			fgets(line, sizeof(line), fp);
			if(!isspace(*line) && (*line)) {
				
				int index;
				double onEdgeDegree, Criticality;
				int _concave, _convex, _saddle, _asSkeleton;
				sscanf_s(line, "%d %lf %lf %d %d %d", &index, &onEdgeDegree, &Criticality, &_concave, &_convex, &_saddle);
				
				sensorPool[index].index = index;
				sensorPool[index].onEdge_Degree = onEdgeDegree;
				sensorPool[index].Criticality = Criticality;
				sensorPool[index].concave = _concave ? true:false;
				sensorPool[index].convex = _convex ? true:false;
				sensorPool[index].saddle = _saddle ? true:false;
				/*sensorPool[index].asSkeleton = _asSkeleton ? true:false;*/
				sensorPool[index].field = this;
			}
		}while(!feof(fp));		
		fclose(fp);
	}
	for(int i = 0; i < nSensors; ++i)
	{

		printf("Processing the %dth sensor.\n", i);

		if(sensorPool[i].onEdge_Degree > 0)
		{
			BoundarySensorSet.push_back(&sensorPool[i]);
		}
	}
}

void Field::recoverPreprocessedTriangle(char *topoFile)
{
	TriangleSet.clear();
	if(topoFile) {
		FILE *fp;
		errno_t err;
		err = fopen_s(&fp, topoFile, "r");
		char line[2048] = {0};
		if (err)  {
			fprintf_s(stderr, "no topo file %s found\n", topoFile);
			exit(0);
		}
		do {
			memset(line, 0, sizeof(line));
			fgets(line, sizeof(line), fp);
			if(!isspace(*line) && (*line)) {

				Triangle T;
				sscanf_s(line, "%d, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n", &T.Tri_ID, 
					&T.t1->location.x, &T.t1->location.y, &T.t1->location.z, 
					&T.t2->location.x, &T.t2->location.y, &T.t2->location.z, 
					&T.t3->location.x, &T.t3->location.y, &T.t3->location.z, &T.Normalx, &T.Normaly, &T.Normalz);

				T.Tri_field = this;
				TriangleSet.push_back(T);
			}
		}while(!feof(fp));		
		fclose(fp);
	}
	for (unsigned int i = 0; i < TriangleSet.size(); ++i)
	{
		if(i == 40)
			printf("wait\n");
		TriangleSet[i].EstablishNeighborhood();
	}
}

void Field::dumpTriangle(char * filename)
{
	if (!filename)
		return;
	FILE *fp;
	errno_t err;
	err = fopen_s(&fp, filename, "a");
	if (!err) 
	{
		for (unsigned int i = 0; i < TriangleSet.size(); ++i)
		{
			Triangle *T = &TriangleSet[i];
			fprintf_s(fp, "%d, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n", T->Tri_ID, 
				T->t1->location.x, T->t1->location.y, T->t1->location.z, 
				T->t2->location.x, T->t2->location.y, T->t2->location.z, 
				T->t3->location.x, T->t3->location.y, T->t3->location.z, T->Normalx, T->Normaly, T->Normalz);
		}
	}
	else 
	{
		fprintf_s(stderr, "open file %s to write failed!\n", filename);
		return;
	}
	fclose(fp);
}
void Field::dumpEdges(char * filename)
{
	if (!filename)
		return;
	FILE *fp;
	errno_t err;
	err = fopen_s(&fp, filename, "w");
	if (!err) 
	{
		for (unsigned int i = 0; i < EdgeSet.size(); ++i)
		{
			LineSegment *L = &EdgeSet[i];

			Sensor* s = L->start;
			Sensor* e = L->end;

			if (s->index < e->index)
			{
				fprintf_s(fp, "%d, %d, %d, %d\n",L->LineID, s->index,e->index,L->usedtimes);
			}
			else
				fprintf_s(fp, "%d, %d, %d, %d\n", L->LineID, e->index,s->index,L->usedtimes);
		}
	}
	else 
	{
		fprintf_s(stderr, "open file %s to write failed!\n", filename);
		return;
	}
	fclose(fp);

}
void Field::dumpTriangleWithID(char * filename)
{
	if (!filename)
		return;
	FILE *fp;
	errno_t err;
	err = fopen_s(&fp, filename, "w");
	if (!err) 
	{
		for (unsigned int i = 0; i < TriangleSet.size(); ++i)
		{
			Triangle *T = &TriangleSet[i];
			fprintf_s(fp, "%d, %d, %d, %d, %lf, %lf, %lf\n", T->Tri_ID, 
				T->t1->index, T->t2->index, T->t3->index,
				T->Normalx, T->Normaly, T->Normalz);
		}
	}
	else 
	{
		fprintf_s(stderr, "open file %s to write failed!\n", filename);
		return;
	}
	fclose(fp);
}
//for debugging circle
void Field::dumpPath(char * filename, Path Depath)
{
	if (!filename)
		return;
	FILE *fp;
	errno_t err;
	err = fopen_s(&fp, filename, "a");
	if (!err) 
	{
		for(unsigned int i = 0; i < Depath.path.size(); ++i)
		{
			Sensor *dump = Depath.path[i];
			fprintf_s(fp, "%d %lf %lf %d %d %d %d %lf %lf %lf\n", dump->index, dump->onEdge_Degree, dump->Criticality, dump->landmark, dump->convex, dump->concave, dump->VoronoiCellID.size(), dump->location.x, dump->location.y, dump->location.z);
		}
		fprintf_s(fp, "#\n");
	}
	else 
	{
		fprintf_s(stderr, "open file %s to write failed!\n", filename);
		return;
	}
	fclose(fp);
}
/*
//for showing path
void Field::dumpPath(char * filename, Path Depath)
{
	if (!filename)
		return;
	FILE *fp;
	errno_t err;
	err = fopen_s(&fp, filename, "a");
	if (!err) 
	{
		for(unsigned int i = 0; i < Depath.path.size(); ++i)
		{
			Sensor *dump = Depath.path[i];
			if(dump->index == Depath.start->index || dump->index == Depath.end->index)
				fprintf_s(fp, "%d %f %f %d %d %d %d %lf %lf %lf %d\n", dump->index, dump->onEdge_Degree, dump->Criticality, dump->landmark, dump->convex, dump->concave, dump->VoronoiCellID.size(), dump->location.x, dump->location.y, dump->location.z, true);
			else
				fprintf_s(fp, "%d %f %f %d %d %d %d %lf %lf %lf %d\n", dump->index, dump->onEdge_Degree, dump->Criticality, dump->landmark, dump->convex, dump->concave, dump->VoronoiCellID.size(), dump->location.x, dump->location.y, dump->location.z, false);
		}
		
	}
	else 
	{
		fprintf_s(stderr, "open file %s to write failed!\n", filename);
		return;
	}
	fclose(fp);
}
*/

void Field::dumpCircle(char * filename, CircleArea ca)
{
	if (!filename)
		return;
	FILE *fp;
	errno_t err;
	err = fopen_s(&fp, filename, "a");
	if (!err) 
	{
		for(SCITER ASSpos = ca.AreaSensorSet.begin(); ASSpos != ca.AreaSensorSet.end(); ++ASSpos)
		{
			Sensor *dump = *ASSpos;
			fprintf_s(fp, "%d %d %lf %lf %lf %lf %d %d %lf\n", dump->index, dump->on_boundary_circle, dump->Criticality, dump->location.x, dump->location.y, dump->location.z, dump->asCenter, dump->concave, ca.AvgHopdist);
		}

	}
	else 
	{
		fprintf_s(stderr, "open file %s to write failed!\n", filename);
		return;
	}
	fprintf_s(fp, "#\n");
	fclose(fp);
}

void Field::dumpLine(char *filename, Sensor *start, Sensor *end)
{
	if (!filename)
		return;
	FILE *fp;
	errno_t err;
	err = fopen_s(&fp, filename, "a");
	if (!err) 
	{
		fprintf_s(fp, "%lf, %lf, %lf, %lf, %lf, %lf\n", 
			start->location.x, start->location.y, start->location.z, 
			end->location.x, end->location.y, end->location.z);
	}
	else 
	{
		fprintf_s(stderr, "open file %s to write failed!\n", filename);
		return;
	}
	fclose(fp);
}

void Field::ResetCancelState()
{
	for(int i = 0; i < nSensors; ++i)
	{
		if(sensorPool[i].Canceled)
		{
			sensorPool[i].Canceled =  false;
		}
	}
}

bool Field::EndGenerating()
{
	unsigned int i;
	long processednum = 0;
	printf("Total edge number is %d\n", EdgeSet.size());
	if(EdgeSet.size() == 352)
		printf("wait\n");
	for(i = 0; i < EdgeSet.size(); ++i)
	{
		if(EdgeSet[i].usedtimes == 2)
			++processednum;
	}
	RemainingEdges = EdgeSet.size()-processednum;
	printf("Need to processing edges' number: %d\n", RemainingEdges);

	for(i = 0; i < EdgeSet.size(); ++i)
	{
		if(EdgeSet[i].usedtimes < 2 && !EdgeSet[i].ProcessingLater)
		{
			return false;
		}
	}
	return true;
}

bool Field::EndUpdating()
{
	unsigned int i;

	for(i = 0; i < RefineEdgeSet.size(); ++i)
	{
		if(RefineEdgeSet[i]->usedtimes == 2)
		{
			return false;
		}
	}
	return true;
}

#if 0
void Field::ReportTopologyState()
{
	int wrongEdgeNum = 0;
	int needtoprocess = 0;
	for(unsigned int i  = 0; i < EdgeSet.size(); ++i)
	{
		if(EdgeSet[i].usedtimes > 2)
			++wrongEdgeNum;
		if(EdgeSet[i].usedtimes < 2)
			++needtoprocess;
	}
	printf("The size of edge set is %d\n", EdgeSet.size());
	printf("The size of refining edge set is %d\n", RefineEdgeSet.size());
	printf("The number of wrong edge is %d\n", wrongEdgeNum);

	//Computing distance error
	//set<Sensor*> OuterBoundaryNodes;
	double sum_dist_error = 0;
	for(int i = 0; i < nSensors; ++i)
	{
		if(!sensorPool[i].Be_used_in_triangulation /*&& sensorPool[i].onEdge_Degree != 0*/)
		{
			Sensor *s = &sensorPool[i];
			double sum_distance = 0xffff;
			Sensor *t1 = NULL;
			Sensor *t2 = NULL;
			Sensor *t3 = NULL;
			for(unsigned int j = 0; j < TriangleSet.size(); ++j)
			{
				double sd = 0;
				sd += s->location.distance(TriangleSet[j].t1->location);
				sd += s->location.distance(TriangleSet[j].t2->location);
				sd += s->location.distance(TriangleSet[j].t3->location);
				if(sd < sum_distance)
				{
					sum_distance = sd;
					t1 = TriangleSet[j].t1;
					t2 = TriangleSet[j].t2;
					t3 = TriangleSet[j].t3;
				}
			}
			//Get normal of the triangle
			double N1x = t1->location.x - t2->location.x;
			double N1y = t1->location.y - t2->location.y;
			double N1z = t1->location.z - t2->location.z;
			double N2x = t1->location.x - t3->location.x;
			double N2y = t1->location.y - t3->location.y;
			double N2z = t1->location.z - t3->location.z;
			//Normal Vector
			double Nx = N1y/N1x*(N1x*N2z - N2x*N1z)/(N1x*N2y - N2x*N1y) - N1z/N1x;
			double Ny = (N2x*N1z - N1x*N2z)/(N1x*N2y  - N2x*N1y);
			double Nz = 1;
			if (N1z == 0 || N1x == 0 || N1y == 0 || N2z == 0 || N2x == 0 || N2y == 0 || (N1x*N2y-N2x*N1y)==0)
			{
				fprintf_s(stdout, "Vector problem!\n");
				exit(0);
			}
			//confirm the notation of normal
			int xrange = 200;
			int yrange = 200;
			int zrange = 200;
			Point Orientation_center;
			
			//for S topology
			if(s->location.y>=-5 && s->location.y <= yrange/10*3 && s->location.z >=-5 && s->location.z<=zrange+5 && s->location.x>=-5 && s->location.x<=xrange+5)
			{
				Orientation_center.x = xrange/2;
				Orientation_center.y = yrange/20*3;
				Orientation_center.z = zrange/2;
			}
			else if(s->location.y>=yrange/10*3 && s->location.y<=yrange/10*7 && s->location.z>=-5 && s->location.z<=zrange+5 && s->location.x>=-5 && s->location.x<=xrange+5)
			{
				Orientation_center.x = xrange/2;
				Orientation_center.y = yrange/2;
				Orientation_center.z = zrange/2;
			}
			else if(s->location.y>=yrange/10*7 && s->location.y<=yrange+5 && s->location.z>=-5 && s->location.z<=zrange+5 && s->location.x>=-5 && s->location.x<=xrange+5)
			{
				Orientation_center.x = xrange/2;
				Orientation_center.y = yrange/20*17;
				Orientation_center.z = zrange/2;
			}
			else
			{
				Orientation_center.x = 0;
				Orientation_center.y = 0;
				Orientation_center.z = 0;
			}

			//for H topology
			//for Y topology
			double reference_vector_x = t1->location.x-Orientation_center.x;
			double reference_vector_y = t1->location.y-Orientation_center.y;
			double reference_vector_z = t1->location.z-Orientation_center.z;
			double product_value = reference_vector_x*Nx + reference_vector_y*Ny + reference_vector_z*Nz;
			if(product_value < 0)
			{
				Nx = 0-Nx;
				Ny = 0-Ny;
				Nz = 0-Nz;
			}
			double ve_x = s->location.x-t1->location.x;
			double ve_y = s->location.y-t1->location.y;
			double ve_z = s->location.z-t1->location.z;
			double ve_product = ve_x*Nx + ve_y*Ny + ve_z*Nz;
			if(ve_product > 0)
			{
				s->Out = true;
				++num_of_outer_boundary_surface;
				double A = fabs(Nx*(s->location.x-t1->location.x) + Ny*(s->location.y - t1->location.y) + Nz*(s->location.z - t1->location.z));
				double B = sqrt(Nx*Nx+Ny*Ny+Nz*Nz);
				double d = A/B;
				sum_dist_error += d;
			}
		}
	}
	Avg_distance_err = sum_dist_error/num_of_outer_boundary_surface;
}
#endif

void Field::ReportTopologyState()
{
	int wrongEdgeNum = 0;
	int needtoprocess = 0;
	for(unsigned int i  = 0; i < EdgeSet.size(); ++i)
	{
		if(EdgeSet[i].usedtimes > 2)
			++wrongEdgeNum;
		if(EdgeSet[i].usedtimes < 2)
			++needtoprocess;
	}
	printf("The size of edge set is %d\n", EdgeSet.size());
	printf("The size of refining edge set is %d\n", RefineEdgeSet.size());
	printf("The number of wrong edge is %d\n", wrongEdgeNum);

	
	//Computing distance error
	//set<Sensor*> OuterBoundaryNodes;
	double sum_dist_error = 0;
	//refine triangles' normal vector
	for(unsigned int i = 0; i < TriangleSet.size(); ++i)
	{
		Triangle *t = &TriangleSet[i];
		Point reference_point;
		reference_point.x = (t->t1->location.x+t->t2->location.x+t->t3->location.x)/3;
		reference_point.y = (t->t1->location.y+t->t2->location.y+t->t3->location.y)/3;
		reference_point.z = (t->t1->location.z+t->t2->location.z+t->t3->location.z)/3;

		//confirm the notation of normal
		int xrange = 200;
		int yrange = 200;
		int zrange = 200;
		Point Orientation_center;
		int Length = LENGTH;
		Orientation_center.x = Length*t->Normalx+reference_point.x;
		Orientation_center.y = Length*t->Normaly+reference_point.y;
		Orientation_center.z = Length*t->Normalz+reference_point.z;
		/*
		//for S topology
		if((Orientation_center.y>=0 && Orientation_center.y <= yrange/5 && Orientation_center.z>=0 && Orientation_center.z<=zrange && Orientation_center.x>=0 && Orientation_center.x<=xrange) 
			|| (Orientation_center.y>=yrange/5*2 && Orientation_center.y<=yrange/5*3 && Orientation_center.z>=0 && Orientation_center.z<=zrange && Orientation_center.x>=0 && Orientation_center.x<=xrange) 
			|| (Orientation_center.y>=yrange/5*4 && Orientation_center.y<=yrange && Orientation_center.z>=0 && Orientation_center.z<=zrange && Orientation_center.x>=0 && Orientation_center.x<=xrange)
			|| (Orientation_center.y>=yrange/5 && Orientation_center.y<=yrange/5*2 && Orientation_center.x>=xrange/5*3 && Orientation_center.x<=xrange && Orientation_center.z>=0 && Orientation_center.z<=zrange) 
			|| (Orientation_center.x>=0 && Orientation_center.x<=xrange/5*2 && Orientation_center.y>=yrange/5*3 && Orientation_center.y<=yrange/5*4 && Orientation_center.z>=0 && Orientation_center.z<=zrange))
		{
			t->Normalx = 0-t->Normalx;
			t->Normaly = 0-t->Normaly;
			t->Normalz = 0-t->Normalz;
		}
		*/
		//for chicago airport & cube with H
		if(WithinTopology(Orientation_center))
		{
			t->Normalx = 0-t->Normalx;
			t->Normaly = 0-t->Normaly;
			t->Normalz = 0-t->Normalz;
		}
	}
	//refine again
	for(unsigned int i = 0; i < TriangleSet.size(); ++i)
	{
		Triangle *t = &TriangleSet[i];
		int matchtimes = 0;
		if(t->TriangleNeighbors.size() > 3)
		{
			matchtimes -= (t->TriangleNeighbors.size()-3);
		}
		for(set<Triangle*>::iterator TNpos = t->TriangleNeighbors.begin(); TNpos != t->TriangleNeighbors.end(); ++TNpos)
		{
			Triangle *tn = *TNpos;
			double product = t->Normalx*tn->Normalx + t->Normaly*tn->Normaly + t->Normalz*tn->Normalz;
			if(product < 0)
			{
				++matchtimes;
			}
		}
		if(matchtimes > 2)
		{
			t->Normalx = 0-t->Normalx;
			t->Normaly = 0-t->Normaly;
			t->Normalz = 0-t->Normalz;
		}
	}
	//////////////////////////////////////////////////////////////////////////
	for(int i = 0; i < nSensors; ++i)
	{
		if(!sensorPool[i].Be_used_in_triangulation /*&& sensorPool[i].onEdge_Degree != 0*/)
		{
			Sensor *s = &sensorPool[i];
			double sum_distance = 0xffff;
			Triangle *t = NULL;
			for(unsigned int j = 0; j < TriangleSet.size(); ++j)
			{
				double sd = 0;
				sd += s->location.distance(TriangleSet[j].t1->location);
				sd += s->location.distance(TriangleSet[j].t2->location);
				sd += s->location.distance(TriangleSet[j].t3->location);
				if(sd < sum_distance)
				{
					sum_distance = sd;
					t = &TriangleSet[j];
				}
			}
			

			Point reference_point;
			reference_point.x = (t->t1->location.x+t->t2->location.x+t->t3->location.x)/3;
			reference_point.y = (t->t1->location.y+t->t2->location.y+t->t3->location.y)/3;
			reference_point.z = (t->t1->location.z+t->t2->location.z+t->t3->location.z)/3;

			double ve_x = s->location.x-reference_point.x;
			double ve_y = s->location.y-reference_point.y;
			double ve_z = s->location.z-reference_point.z;
			double ve_product = ve_x*t->Normalx + ve_y*t->Normaly + ve_z*t->Normalz;
			if(ve_product > 0)
			{
				s->Out = true;
				++num_of_outer_boundary_surface;
				double A = fabs(t->Normalx*(s->location.x-t->t1->location.x) + t->Normaly*(s->location.y - t->t1->location.y) + t->Normalz*(s->location.z - t->t1->location.z));
				double B = sqrt(t->Normalx*t->Normalx+t->Normaly*t->Normaly+t->Normalz*t->Normalz);
				double d = A/B;
				sum_dist_error += d;
			}
		}
	}
	Avg_distance_err = sum_dist_error/nSensors;
}

void Field::DetectEdgeStateforSkeleton()
{
	//Establish the edge set
	int LineID = 0;
	for (unsigned int i = 0; i < LandmarkSet.size(); ++i)
	{
		Sensor *l = LandmarkSet[i];
		for(set<int>::iterator lcdpos = l->ConnectedID.begin(); lcdpos != l->ConnectedID.end(); ++lcdpos)
		{
			int lcID = *lcdpos;
			if(EdgeSet.empty())
			{
				LineSegment edge(0xffff, l, &sensorPool[lcID]);
				edge.LineID = LineID;
				++LineID;
				EdgeSet.push_back(edge);
			}
			else
			{
				LineSegment edge(0xffff, l, &sensorPool[lcID]);
				bool sameline = false;
				for(unsigned int j = 0; j < EdgeSet.size(); ++j)
				{
					if(edge.TheSameLine(&EdgeSet[j]))
					{
						sameline = true;
						break;
					}
				}
				if(!sameline)
				{
					edge.LineID = LineID;
					++LineID;
					EdgeSet.push_back(edge);
				}
			}
		}
	}
	//detect triangles
	for(unsigned int i = 0; i < EdgeSet.size(); ++i)
	{
		set<Sensor*> candidate3rdpointsSet;
		set<Sensor*> filterSet;
		for(unsigned int j = 0; j < LandmarkSet.size(); ++j)
		{
			Sensor *ls = LandmarkSet[j];
			int match_times = 0;
			if(!EdgeSet[i].CanFindASensorInThe3rdPoint(ls))
			{
				for(set<int>::iterator lcdpos = ls->ConnectedID.begin(); lcdpos != ls->ConnectedID.end(); ++lcdpos)
				{
					int lcID = *lcdpos;
					if(lcID == EdgeSet[i].start->index || lcID == EdgeSet[i].end->index)
					{
						match_times++;
					}
				}
				if(match_times == 2)
				{
					candidate3rdpointsSet.insert(ls);
					//dumpTriangle("triangle", EdgeSet[i].start, EdgeSet[i].end, ls);
					//TriangleSet.push_back(Triangle(EdgeSet[i].start, EdgeSet[i].end, ls));
				}
			}
		}
		if(candidate3rdpointsSet.size() == 2)
		{
			for(SCITER cpspos = candidate3rdpointsSet.begin(); cpspos != candidate3rdpointsSet.end(); ++cpspos)
			{
				Sensor *s = *cpspos;
				EdgeSet[i].The3rdPoint.push_back(s);
				EdgeSet[i].usedtimes++;
			}
		}
		else
		{
			double distance =  0xffff;
			Sensor *ChooseSensor = NULL;
			for(SCITER cpspos1 = candidate3rdpointsSet.begin(); cpspos1 != candidate3rdpointsSet.end(); ++cpspos1)
			{
				bool connected = false;
				Sensor *s1 = *cpspos1;
				for (set<int>::iterator CIDpos = s1->ConnectedID.begin(); CIDpos != s1->ConnectedID.end(); ++CIDpos)
				{
					int CID = *CIDpos;
					for(SCITER cpspos2 = candidate3rdpointsSet.begin(); cpspos2 != candidate3rdpointsSet.end(); ++cpspos2)
					{
						Sensor *s2 = *cpspos2;
						if(s1 != s2)
						{
							if(CID == s2->index)
							{
								filterSet.insert(s1);
								filterSet.insert(s2);
								for(SCITER fspos = filterSet.begin(); fspos != filterSet.end(); ++fspos)
								{
									Sensor *fs = *fspos;
									double distance0 = fs->location.distance(EdgeSet[i].start->location)+fs->location.distance(EdgeSet[i].end->location);
									if(distance > distance0)
									{
										distance = distance0;
										ChooseSensor = fs;
									}
								}
								if(!EdgeSet[i].CanFindASensorInThe3rdPoint(ChooseSensor))
								{
									EdgeSet[i].The3rdPoint.push_back(ChooseSensor);
									EdgeSet[i].usedtimes++;
								}
								filterSet.clear();
								connected = true;
							}
						}
					}
				}
				if (!connected)
				{
					EdgeSet[i].The3rdPoint.push_back(s1);
					EdgeSet[i].usedtimes++;
				}
			}
		}
	}
	//Edge flip
	vector<LineSegment*> EdgeflipSet;
	for(unsigned int i = 0; i < EdgeSet.size(); ++i)
	{
		if(EdgeSet[i].usedtimes > 2)
		{
			EdgeflipSet.push_back(&EdgeSet[i]);
		}
		if(EdgeSet[i].usedtimes < 2)
		{
			printf("Warning!!!Edge %d used only once\n");
		}
		if(EdgeSet[i].usedtimes == 2)
		{
			FinalEdgeSet.push_back(EdgeSet[i]);
		}
	}
	for(unsigned int i = 0; i < EdgeflipSet.size(); ++i)
	{
		int dist = 0;
		LineSegment Forbidden;
		for(unsigned int j = 0; j < EdgeflipSet[i]->The3rdPoint.size(); ++j)
		{
			Sensor *start = EdgeflipSet[i]->The3rdPoint[j];
			for(unsigned int k = j; k < EdgeflipSet[i]->The3rdPoint.size(); ++k)
			{
				Sensor *end = EdgeflipSet[i]->The3rdPoint[k];
				if(start != end)
				{
					int hop = start->FloodingforGetHopsBetweenStartandEndforSkeleton(end);
					if(dist < hop)
					{
						Forbidden.start = start;
						Forbidden.end = end;
					}
				}
			}
		}
		for(unsigned int j = 0; j < EdgeflipSet[i]->The3rdPoint.size(); ++j)
		{
			Sensor *start = EdgeflipSet[i]->The3rdPoint[j];
			for(unsigned int k = j; k < EdgeflipSet[i]->The3rdPoint.size(); ++k)
			{
				Sensor *end = EdgeflipSet[i]->The3rdPoint[k];
				if(start != end)
				{
					LineSegment l(0xffff, start, end);
					if(!l.TheSameLine(&Forbidden))
					{
						l.The3rdPoint.push_back(EdgeflipSet[i]->start);
						l.The3rdPoint.push_back(EdgeflipSet[i]->end);
						l.usedtimes = l.The3rdPoint.size();
						FinalEdgeSet.push_back(l);
					}
				}
			}
		}
	}
	//distribute line ID & dump triangles
	int Tri_ID = 0;
	for(unsigned int i = 0; i < FinalEdgeSet.size(); ++i)
	{
		FinalEdgeSet[i].LineID = i;
		for(unsigned int j = 0; j < FinalEdgeSet[i].The3rdPoint.size(); ++j)
		{
			Triangle T(Tri_ID, FinalEdgeSet[i].start, FinalEdgeSet[i].end, FinalEdgeSet[i].The3rdPoint[j]);
			T.Tri_field = this;
			//dumpTriangle("triangle");
			++Tri_ID;
		}
	}
	int wrongEdgeNum = 0;
	int needtoprocess = 0;
	for(unsigned int i  = 0; i < EdgeSet.size(); ++i)
	{
		if(EdgeSet[i].usedtimes > 2)
			++wrongEdgeNum;
		if(EdgeSet[i].usedtimes < 2)
			++needtoprocess;
	}
	printf("The size of edge set is %d\n", EdgeSet.size());
	printf("The size of refining edge set is %d\n", needtoprocess);
	printf("The number of wrong edge is %d\n", wrongEdgeNum);
}

bool Field::LandmarkBelongToRefineSensorSet(Sensor *l)
{
	for(SCITER RefineSensorSetpos = RefineSensorSet.begin(); RefineSensorSetpos != RefineSensorSet.end(); ++RefineSensorSetpos)
	{
		Sensor *s = *RefineSensorSetpos;
		if(l->index == s->index)
		{
			return true;
		}
	}
	return false;
}

void Field::Triangulation()
{
	int TriID = 0;
	int debugnum = 0;
	int LineIDindex = 0;
	for(int iter = 0; iter < NumOfBoundaryPart; ++iter)
	{
		int SI = EdgeSet.size();
		int i = 0;
		int j = 0;
		//Triangulation
		bool CancelThisPoint = false;
		bool l1samelinetag = false;
		bool l2samelinetag = false;
		Sensor *start = NULL;
		Sensor *nei = NULL;
		SCITER LMNeipos;
		multimap< double, Sensor*, greater<double> > CommonPoints;
		multimap< double, Sensor*, greater<double> >::iterator CommonPointspos;
		multimap< double, Sensor*, greater<double> > SharePoints;
		for(int la = 0; la < LandmarkSet.size(); ++la)
		{
			Sensor *s = LandmarkSet[la];
			if(!s->picked)
				start = s;
		}
		int current_partID = start->boundary_partID;
		for(int la = 0; la < LandmarkSet.size(); ++la)
		{
			Sensor *s = LandmarkSet[la];
			if(s->boundary_partID == current_partID)
				s->picked = true;
		}
		LMNeipos = start->LandmarkNeighbors.begin();
		nei = *LMNeipos;
		EdgeSet.push_back(LineSegment(LineIDindex, start, nei));
		start->Be_used_in_triangulation = true;
		nei->Be_used_in_triangulation = true;
		EdgeSet[SI].Lpath.EstablishPath(start, nei);
		EdgeSet[SI].Lpath.pafield = this;  
		for(unsigned int p = 1; p < EdgeSet[0].Lpath.path.size()-1; ++p)
		{
			Sensor *lp = EdgeSet[0].Lpath.path[p];
			lp->line_ends_ID.insert(make_pair(start->index, nei->index));
		}
		//EdgeSet[0].start->LandmarkVicilityLines.insert(&EdgeSet[0]);
		//EdgeSet[0].end->LandmarkVicilityLines.insert(&EdgeSet[0]);
		//ConnectingLandmarkIDPairs.insert(make_pair(start->index, nei->index));
		++LineIDindex;
		i = SI;
		int debug_count = 0;
		do
		{
			SCITER LandmarkNeighborspos;
			if(EdgeSet[i].usedtimes > 2)
			{
				printf("Edge used times wrong!\n");
				exit(0);
			}
			while(EdgeSet[i].usedtimes < 2)
			{
				++EdgeSet[i].Processingtimes;
				if(EdgeSet[i].Processingtimes >= 10)
				{
					EdgeSet[i].ProcessingLater = true;
					break;
				}
				Sensor *Greatest3rdPoint = NULL;
				CommonPoints.clear();
				SharePoints.clear();
				//CannotFindASensorInLandmarkNeighbors = true;
				//CrossingOtherLines = true;
				//DismissThe3rdPoint = true;
				for(LandmarkNeighborspos = EdgeSet[i].start->LandmarkNeighbors.begin(); LandmarkNeighborspos != EdgeSet[i].start->LandmarkNeighbors.end(); ++LandmarkNeighborspos)
				{
					Sensor *S0 = *LandmarkNeighborspos;
					if(S0->index != EdgeSet[i].end->index/* The 3rd point is not end of line*/ 
						&& EdgeSet[i].end->CanFindASensorInLandmarkNeighbors(S0) 
						&& !EdgeSet[i].CanFindASensorInThe3rdPoint(S0) 
						&& !EdgeSet[i].CorrespThe3rdPointOfLineVertex(S0) 
						&& !S0->Canceled)
					{
						CommonPoints.insert(make_pair(S0->onEdge_Degree, S0));
					}
				}
				if(!CommonPoints.empty())
				{
					for(CommonPointspos = CommonPoints.begin(); CommonPointspos != CommonPoints.end(); ++CommonPointspos)
					{
						Sensor *S0 = CommonPointspos->second;
						if(noCrossingBy(&EdgeSet[i], S0) 
							&& !S0->DismissThisSensorByHop(&EdgeSet[i], 2)
							&& !S0->DismissThisSensor(TriangleSet))
						{
							SharePoints.insert(make_pair(S0->Rou_p, S0));
						}
					}
				}
				else
				{
					EdgeSet[i].ProcessingLater = true;
					break;
				}
				/*
				else
				{
					for (CommonPointspos = CommonPoints.begin(); CommonPointspos != CommonPoints.end(); ++CommonPointspos)
					{
						Sensor *S0 = CommonPointspos->second;
						SharePoints.insert(make_pair(S0->onEdge_Degree, S0));
					}
				}
				*/
				if(SharePoints.empty())
				{
					for(CommonPointspos = CommonPoints.begin(); CommonPointspos != CommonPoints.end(); ++CommonPointspos)
					{
						Sensor *S0 = CommonPointspos->second;
						if(noCrossingBy(&EdgeSet[i], S0) 
							&& !S0->DismissThisSensorByHop(&EdgeSet[i], 2))
						{
							SharePoints.insert(make_pair(S0->Rou_p, S0));
						}
					}
				}
				
				if(SharePoints.empty())
				{
					EdgeSet[i].ProcessingLater = true;
					break;
				}
				else
				{
					CancelThisPoint = false;
					Greatest3rdPoint = EdgeSet[i].GetTheGreatest3rdPoint(SharePoints);
					if(Greatest3rdPoint == NULL)
					{
						EdgeSet[i].ProcessingLater = true;
					}
					else
					{
						LineSegment *tempLine = NULL;
						l1samelinetag = false;
						LineSegment l1(LineIDindex, EdgeSet[i].end, Greatest3rdPoint);
						
						for(unsigned int j = 0; j < EdgeSet.size(); ++j)
						{
							if(l1.TheSameLine(&EdgeSet[j]))
							{
								l1samelinetag = true;
								if(EdgeSet[j].usedtimes < 2)
								{
									tempLine = &EdgeSet[j];
									break;
								}
								else
								{
									CancelThisPoint = true;
									break;
								}
							}	
						}
						if(!l1samelinetag)
						{
							l1.The3rdPoint.push_back(EdgeSet[i].start);
							l1.usedtimes = l1.The3rdPoint.size();
						}
						l2samelinetag = false;
						LineSegment l2(LineIDindex, Greatest3rdPoint, EdgeSet[i].start);
						
						for(unsigned int j = 0; j < EdgeSet.size(); ++j)
						{
							if(l2.TheSameLine(&EdgeSet[j]))
							{
								l2samelinetag = true;
								if(EdgeSet[j].usedtimes < 2  && !CancelThisPoint)
								{
									EdgeSet[j].The3rdPoint.push_back(EdgeSet[i].end);
									EdgeSet[j].usedtimes = EdgeSet[j].The3rdPoint.size();
									EdgeSet[j].ProcessingThe3rdPoint();
									if(tempLine != NULL)
									{
										tempLine->The3rdPoint.push_back(EdgeSet[i].start);
										tempLine->usedtimes = tempLine->The3rdPoint.size();
										tempLine->ProcessingThe3rdPoint();
									}
									break;
								}
								else
								{
									CancelThisPoint = true;
									break;
								}
							}
						}
						if(!l2samelinetag)
						{
							l2.The3rdPoint.push_back(EdgeSet[i].end);
							l2.usedtimes = l2.The3rdPoint.size();
							if(tempLine != NULL)
							{
								tempLine->The3rdPoint.push_back(EdgeSet[i].start);
								tempLine->usedtimes = tempLine->The3rdPoint.size();
								tempLine->ProcessingThe3rdPoint();
							}
						}
						if(!CancelThisPoint)
						{
							if(!l1samelinetag)
							{
								l1.LineID = LineIDindex;
								l1.Lpath.EstablishPath(l1.start, l1.end);
								l1.Lpath.pafield = this;
								for(unsigned int p = 1; p < l1.Lpath.path.size()-1; ++p)
								{
									Sensor *lp = l1.Lpath.path[p];
									lp->line_ends_ID.insert(make_pair(l1.start->index, l1.end->index));
								}
								EdgeSet.push_back(l1);
								
								++LineIDindex;
							}
							if(!l2samelinetag)
							{
								l2.LineID = LineIDindex;
								l2.Lpath.EstablishPath(l2.start, l2.end);
								l2.Lpath.pafield = this;
								for(unsigned int p = 1; p < l2.Lpath.path.size()-1; ++p)
								{
									Sensor *lp = l2.Lpath.path[p];
									lp->line_ends_ID.insert(make_pair(l2.start->index, l2.end->index));
								}
								EdgeSet.push_back(l2);
								
								++LineIDindex;
							}
							EdgeSet[i].The3rdPoint.push_back(Greatest3rdPoint);
							EdgeSet[i].usedtimes = EdgeSet[i].The3rdPoint.size();
							EdgeSet[i].ProcessingThe3rdPoint();
							
							EdgeSet[i].start->Be_used_in_triangulation = true;
							EdgeSet[i].end->Be_used_in_triangulation = true;
							Greatest3rdPoint->Be_used_in_triangulation = true;
							Triangle T(TriID, EdgeSet[i].start, EdgeSet[i].end, Greatest3rdPoint);
							T.Tri_field = this;
							T.GetNormalVector();
							TriangleSet.push_back(T);
							++debug_count;
							//dumpTriangle("triangle", &T);
							++TriID;
							if(debug_count == 59)
								bool stop = 1;
						}
						else
						{
							//printf("Processing the %dth line, Cenceled the sensor %d\n", i, Greatest3rdPoint->index );
							Greatest3rdPoint->Canceled = true;
							printf("Cancel this point, process again\n");
						}
					}
				}
			}
			++i;
			ResetCancelState();
		} while(!EndGenerating());
		cout << "Iteration " << iter << " complete" << endl;
	}
	
	printf("First step end!\n");
}
#if 0
void Field::TriangulationForSkeleton()
{
	int i = 0;
	int j = 0;

	int debugnum = 0;

	int LineIDindex = 0;

	//Triangulation
	bool CancelThisPoint = false;
	bool l1samelinetag = false;
	bool l2samelinetag = false;
	Sensor *start = NULL;
	Sensor *nei = NULL;
	SCITER LMNeipos;
	multimap< double, Sensor*, greater<double> > CommonPoints;
	multimap< double, Sensor*, greater<double> >::iterator CommonPointspos;
	multimap< double, Sensor*, greater<double> > SharePoints;
	start = LandmarkSet[0];
	LMNeipos = LandmarkSet[0]->LandmarkNeighbors.begin();
	nei = *LMNeipos;
	EdgeSet.push_back(LineSegment(LineIDindex, start, nei));
	start->Be_used_in_triangulation = true;
	nei->Be_used_in_triangulation = true;

	EdgeSet[0].Lpath.EstablishDirectPath(start, nei);
	EdgeSet[0].Lpath.pafield = this;  
	for(unsigned int p = 1; p < EdgeSet[0].Lpath.path.size()-1; ++p)
	{
		Sensor *lp = EdgeSet[0].Lpath.path[p];
		lp->line_ends_ID.insert(make_pair(start->index, nei->index));
	}
	//EdgeSet[0].start->LandmarkVicilityLines.insert(&EdgeSet[0]);
	//EdgeSet[0].end->LandmarkVicilityLines.insert(&EdgeSet[0]);

	//ConnectingLandmarkIDPairs.insert(make_pair(start->index, nei->index));
	++LineIDindex;
	i = 0;
	do
	{

		SCITER LandmarkNeighborspos;

		if(EdgeSet[i].usedtimes > 2)
		{
			printf("Edge used times wrong!\n");
			exit(0);
		}

		while(EdgeSet[i].usedtimes < 2)
		{
			++EdgeSet[i].Processingtimes;
			if(EdgeSet[i].Processingtimes >= 10)
			{
				EdgeSet[i].ProcessingLater = true;

				break;
			}
			Sensor *Greatest3rdPoint = NULL;
			if(i==258)
				printf("wait\n");

			CommonPoints.clear();
			SharePoints.clear();

			//CannotFindASensorInLandmarkNeighbors = true;
			//CrossingOtherLines = true;
			//DismissThe3rdPoint = true;

			for(LandmarkNeighborspos = EdgeSet[i].start->LandmarkNeighbors.begin(); LandmarkNeighborspos != EdgeSet[i].start->LandmarkNeighbors.end(); ++LandmarkNeighborspos)
			{
				Sensor *S0 = *LandmarkNeighborspos;
				if(S0->index != EdgeSet[i].end->index/* The 3rd point is not end of line*/ 
					&& EdgeSet[i].end->CanFindASensorInLandmarkNeighbors(S0) 
					&& !EdgeSet[i].CanFindASensorInThe3rdPoint(S0) 
					&& !EdgeSet[i].CorrespThe3rdPointOfLineVertex(S0) 
					&& !S0->Canceled)
				{
					CommonPoints.insert(make_pair(S0->onEdge_Degree, S0));
				}

			}



			if(!CommonPoints.empty())
			{
				for(CommonPointspos = CommonPoints.begin(); CommonPointspos != CommonPoints.end(); ++CommonPointspos)
				{
					Sensor *S0 = CommonPointspos->second;
					if(noCrossingByforSkeleton(&EdgeSet[i], S0) 
						&& !S0->DismissThisSensorByHop(&EdgeSet[i], 2)
						&& !S0->DismissThisSensor(TriangleSet))
					{
						SharePoints.insert(make_pair(S0->Rou_p, S0));
					}
				}

			}
			else
			{
				EdgeSet[i].ProcessingLater = true;
				break;
			}
			/*
			else
			{
			for (CommonPointspos = CommonPoints.begin(); CommonPointspos != CommonPoints.end(); ++CommonPointspos)
			{
			Sensor *S0 = CommonPointspos->second;
			SharePoints.insert(make_pair(S0->onEdge_Degree, S0));
			}
			}
			*/


			if(SharePoints.empty())
			{
				for(CommonPointspos = CommonPoints.begin(); CommonPointspos != CommonPoints.end(); ++CommonPointspos)
				{
					Sensor *S0 = CommonPointspos->second;
					if(noCrossingByforSkeleton(&EdgeSet[i], S0) 
						&& !S0->DismissThisSensorByHop(&EdgeSet[i], 2))
					{
						SharePoints.insert(make_pair(S0->Rou_p, S0));
					}
				}
			}

			if(SharePoints.empty())
			{
				EdgeSet[i].ProcessingLater = true;
				break;
			}
			else
			{
				CancelThisPoint = false;
				Greatest3rdPoint = EdgeSet[i].GetTheGreatest3rdPointforSkeleton(SharePoints);
				if(Greatest3rdPoint == NULL)
				{
					EdgeSet[i].ProcessingLater = true;
				}
				else
				{
					LineSegment *tempLine = NULL;
					l1samelinetag = false;
					LineSegment l1(LineIDindex, EdgeSet[i].end, Greatest3rdPoint);

					for(unsigned int j = 0; j < EdgeSet.size(); ++j)
					{
						if(l1.TheSameLine(&EdgeSet[j]))
						{
							l1samelinetag = true;
							if(EdgeSet[j].usedtimes < 2)
							{
								tempLine = &EdgeSet[j];
								break;
							}
							else
							{
								CancelThisPoint = true;
								break;
							}
						}	
					}
					if(!l1samelinetag)
					{
						l1.The3rdPoint.push_back(EdgeSet[i].start);
						l1.usedtimes = l1.The3rdPoint.size();
					}
					l2samelinetag = false;
					LineSegment l2(LineIDindex, Greatest3rdPoint, EdgeSet[i].start);

					for(unsigned int j = 0; j < EdgeSet.size(); ++j)
					{
						if(l2.TheSameLine(&EdgeSet[j]))
						{
							l2samelinetag = true;
							if(EdgeSet[j].usedtimes < 2  && !CancelThisPoint)
							{
								EdgeSet[j].The3rdPoint.push_back(EdgeSet[i].end);
								EdgeSet[j].usedtimes = EdgeSet[j].The3rdPoint.size();
								EdgeSet[j].ProcessingThe3rdPoint();
								if(tempLine != NULL)
								{
									tempLine->The3rdPoint.push_back(EdgeSet[i].start);
									tempLine->usedtimes = tempLine->The3rdPoint.size();
									tempLine->ProcessingThe3rdPoint();
								}
								break;
							}
							else
							{
								CancelThisPoint = true;
								break;
							}
						}
					}
					if(!l2samelinetag)
					{
						l2.The3rdPoint.push_back(EdgeSet[i].end);
						l2.usedtimes = l2.The3rdPoint.size();
						if(tempLine != NULL)
						{
							tempLine->The3rdPoint.push_back(EdgeSet[i].start);
							tempLine->usedtimes = tempLine->The3rdPoint.size();
							tempLine->ProcessingThe3rdPoint();
						}
					}
					if(!CancelThisPoint)
					{
						if(!l1samelinetag)
						{
							l1.LineID = LineIDindex;
							l1.Lpath.EstablishDirectPath(l1.start, l1.end);
							l1.Lpath.pafield = this;
							for(unsigned int p = 1; p < l1.Lpath.path.size()-1; ++p)
							{
								Sensor *lp = l1.Lpath.path[p];
								lp->line_ends_ID.insert(make_pair(l1.start->index, l1.end->index));
							}
							EdgeSet.push_back(l1);

							++LineIDindex;
						}
						if(!l2samelinetag)
						{
							l2.LineID = LineIDindex;
							l2.Lpath.EstablishDirectPath(l2.start, l2.end);
							l2.Lpath.pafield = this;
							for(unsigned int p = 1; p < l2.Lpath.path.size()-1; ++p)
							{
								Sensor *lp = l2.Lpath.path[p];
								lp->line_ends_ID.insert(make_pair(l2.start->index, l2.end->index));
							}
							EdgeSet.push_back(l2);

							++LineIDindex;
						}
						EdgeSet[i].The3rdPoint.push_back(Greatest3rdPoint);
						EdgeSet[i].usedtimes = EdgeSet[i].The3rdPoint.size();
						EdgeSet[i].ProcessingThe3rdPoint();
						dumpTriangle("triangle", EdgeSet[i].start, EdgeSet[i].end, Greatest3rdPoint);
						EdgeSet[i].start->Be_used_in_triangulation = true;
						EdgeSet[i].end->Be_used_in_triangulation = true;
						Greatest3rdPoint->Be_used_in_triangulation = true;
						++debugnum;
						if(debugnum == 184)
							printf("wait\n");
						TriangleSet.push_back(Triangle(EdgeSet[i].start, EdgeSet[i].end, Greatest3rdPoint));
					}
					else
					{

						//printf("Processing the %dth line, Cenceled the sensor %d\n", i, Greatest3rdPoint->index );
						Greatest3rdPoint->Canceled = true;
						printf("Cancel this point, process again\n");
					}
				}
			}
		}
		++i;
		ResetCancelState();
	} while(!EndGenerating());
	printf("First step end!\n");
}
#endif

void Field::TriangulationForSkeleton()
{
	int i = 0;
	int j = 0;

	int debugnum = 0;

	int LineIDindex = 0;
	

	//Triangulation
	//CDG
	vector<Sensor*>::iterator LandmarkSetpos;
	bool CancelTheLine = false;
	for(LandmarkSetpos = LandmarkSet.begin(); LandmarkSetpos != LandmarkSet.end(); ++LandmarkSetpos)
	{
		Sensor *CurrentLandmark = *LandmarkSetpos;
		for(SCITER LNpos = CurrentLandmark->LandmarkNeighbors.begin(); LNpos != CurrentLandmark->LandmarkNeighbors.end(); ++LNpos)
		{
			vector<int> AssociateID;
			CancelTheLine = false;
			Sensor *LN = *LNpos;
			Path line(CurrentLandmark, LN);
			line.EstablishDirectPath(CurrentLandmark, LN);
			for(unsigned int p = 1; p < line.path.size()-1; ++p)
			{
				Sensor *l = line.path[p];
				if(l->VoronoiCellID.size()==2)
				{
					set<int>::iterator VIDpos = l->VoronoiCellID.begin();
					int a = *VIDpos;
					++VIDpos;
					int b = *VIDpos;
					if((a == CurrentLandmark->index && b == LN->index) || (a == LN->index && b == CurrentLandmark->index))
						continue;
				}

				for(set<int>::iterator VIDpos = l->VoronoiCellID.begin(); VIDpos !=l->VoronoiCellID.end(); ++VIDpos)
				{
					int VID = *VIDpos;
					if(VID != LN->index && VID != CurrentLandmark->index)
					{
						CancelTheLine = true;
						AssociateID.clear();
						break;
					}
					else
					{
						AssociateID.push_back(VID);
					}
				}
				if(CancelTheLine)
				{
					break;
				}
			}
			if(!CancelTheLine)
			{
				CurrentLandmark->ConnectedID.insert(LN->index);
				LN->ConnectedID.insert(CurrentLandmark->index);
				dumpLine("Line", CurrentLandmark, LN);
				for(unsigned int p = 1; p < line.path.size()-1; ++p)
				{
					Sensor *lp = line.path[p];
					lp->line_ends_ID.insert(make_pair(CurrentLandmark->index, LN->index));
				}
			}
		}
	}
	printf("Skeleton triangulation part end!\n");
}

void Field::EstablishLandmarkVicilitylines()
{
	unsigned int i;
	for(i = 0; i < LandmarkSet.size(); ++i)
	{
		LandmarkSet[i]->LandmarkVicilityLines.clear();
	}
	for(i = 0; i < EdgeSet.size(); ++i)
	{
		Sensor *estart = EdgeSet[i].start;
		Sensor *eend = EdgeSet[i].end;
		estart->LandmarkVicilityLines.insert(&EdgeSet[i]);
		eend->LandmarkVicilityLines.insert(&EdgeSet[i]);
	}
}

bool Field::noCrossingBy(LineSegment *baseline, Sensor* s)
{
	if(baseline->The3rdPoint.empty())
	{
		//CrossingOtherLines = false;
		return true;
	}
	if(baseline->The3rdPoint.size() == 1)
	{
		//EstablishLandmarkVicilitylines();
		Path paths1, paths2;
		paths1.EstablishPath(s, baseline->start);					//waiting for judging
		paths1.pafield = this;
		paths2.EstablishPath(s, baseline->end);						//waiting for judging
		paths2.pafield = this;
		Path pathcors1(baseline->The3rdPoint[0], baseline->end);	//existing path
		pathcors1.pafield = this;
		Path pathcors2(baseline->The3rdPoint[0], baseline->start);	//existing path
		pathcors2.pafield = this;
		//if(paths1.PathCrossing(&pathcors1, 1) || paths2.PathCrossing(&pathcors2, 1))

		if( paths1.PathCrossing(1) || paths2.PathCrossing(1) )
		{
			printf("Detect Crossing!!!\n");
			//CrossingOtherLines = true;
			return false;
		}
		else
		{
			//CrossingOtherLines = false;
			return true;
		}
	}
	else
	{
		printf("Processing Edge problem\n");
		exit(0);
	}
}

bool Field::noCrossingByforSkeleton(LineSegment *baseline, Sensor* s)
{
	if(baseline->The3rdPoint.empty())
	{
		//CrossingOtherLines = false;
		return true;
	}
	if(baseline->The3rdPoint.size() == 1)
	{
		//EstablishLandmarkVicilitylines();
		Path paths1, paths2;
		paths1.EstablishDirectPath(s, baseline->start);					//waiting for judging
		paths1.pafield = this;
		paths2.EstablishDirectPath(s, baseline->end);						//waiting for judging
		paths2.pafield = this;
		Path pathcors1(baseline->The3rdPoint[0], baseline->end);	//existing path
		pathcors1.pafield = this;
		Path pathcors2(baseline->The3rdPoint[0], baseline->start);	//existing path
		pathcors2.pafield = this;
		//if(paths1.PathCrossing(&pathcors1, 1) || paths2.PathCrossing(&pathcors2, 1))

		if( paths1.PathCrossing(1) || paths2.PathCrossing(1) )
		{
			printf("Detect Crossing!!!\n");
			//CrossingOtherLines = true;
			return false;
		}
		else
		{
			//CrossingOtherLines = false;
			return true;
		}
	}
	else
	{
		printf("Processing Edge problem\n");
		exit(0);
	}
}

void Field::DetectSkeleton(int khop, int offset)
{
	//distribute the DistanceToLandmarks to every sensor
	for (unsigned int i = 0; i < LandmarkSet.size(); ++i)
	{
		Sensor *l = LandmarkSet[i];
		l->LocalFloodingForMeasureTheHopsToSensors(khop);
	}

	for(int i = 0; i < nSensors; ++i)
	{
		if(sensorPool[i].onEdge_Degree == 0)
		{
			Sensor *s = &sensorPool[i];
			if(s->DistanceToLandmarks.size() > 1)
			{
				multimap<int, int, less<int> >::iterator DTLpos;
				multimap<int, int, less<int> >::iterator tempDTLpos;
				for(DTLpos = s->DistanceToLandmarks.begin(); DTLpos != s->DistanceToLandmarks.end(); ++DTLpos)
				{
					tempDTLpos = ++DTLpos;
					--DTLpos;
					int v1 = DTLpos->first;
					int ID1 = DTLpos->second;
					int v2 = 0;
					int ID2 = 0xffff;
					if(tempDTLpos != s->DistanceToLandmarks.end())
					{
						v2 = tempDTLpos->first;
						ID2 = tempDTLpos->second;
					}
					if (abs(v1 - v2) <= offset && NotLandmarkNeighbors(ID1, ID2) && PathCrossingInnerNode(ID1, ID2))
					{
						s->asSkeleton = true;
						printf("Detect a skelton point!\n");
					}
					break;
				}
			}
		}
	}
	int num_skeleton = 0;
	for (int i = 0; i < nSensors; ++i)
	{
		if (sensorPool[i].asSkeleton)
		{
			SkeletonNodeSet.push_back(&sensorPool[i]);
			++num_skeleton;
		}
	}
	printf("Total skeleton nodes: %d\n", num_skeleton);
}

bool Field::NotLandmarkNeighbors(int id1, int id2)
{
	for(SCITER LNpos = sensorPool[id1].LandmarkNeighbors.begin(); LNpos != sensorPool[id1].LandmarkNeighbors.end(); ++LNpos)
	{
		Sensor *ln = *LNpos;
		if(ln->index == id2)
		{
			return false;
		}
	}
	return true;
}

bool Field::PathCrossingInnerNode(int id1, int id2)
{
	Sensor *s1 = &sensorPool[id1];
	Sensor *s2 = &sensorPool[id2];
	int count = 0;
	Path s;
	s.EstablishDirectPath(s1, s2);
	for (unsigned int i = 0; i < s.path.size(); ++i)
	{
		Sensor *p = s.path[i];
		if(p->onEdge_Degree == 0)
		{
			count++;
		}
	}
	if(count > (double)s.path.size()/3*1.8)
		return true;
	else
		return false;
}

bool Field::WithinTopology(Point Ore)
{
	double min_dist = 0xffff;
	for (unsigned int i = 0; i < nSensors; ++i)
	{
		double dist = sensorPool[i].location.distance(Ore);
		if(dist < min_dist)
		{
			min_dist = dist;
		}
	}
	if(min_dist <= COMM_RANGE)
		return true;
	else
		return false;
}

void Field::PartitionBoundarySamples()
{
	int partID = 0;
	int count_covered = 0;
	for (int i = 0; i < nSensors; ++i)
	{
		sensorPool[i].broadcastCovered = 0;
	}
	for(unsigned int i = 0 ; i < LandmarkSet.size(); ++i)
	{
		Sensor *s = LandmarkSet[i];
		vector<Sensor*> SameIDSet;
		if(!s->broadcastCovered)
		{
			SameIDSet = s->SpreadPartID(partID);
			if(SameIDSet.size() < 10)
			{
				s->landmark = false;
				for(unsigned int j = 0; j < SameIDSet.size(); ++j)
				{
					Sensor *ss = SameIDSet[j];
					ss->landmark = false;
				}
			}
			else
			{
				++count_covered;
				++partID;
			}
		}
		else
			++count_covered;
	}
	//++partID;
	NumOfBoundaryPart = partID;
	LandmarkSet.clear();
	for(int i = 0; i < nSensors; ++i)
	{
		if(sensorPool[i].landmark)
			LandmarkSet.push_back(&sensorPool[i]);
	}
	cout << "Spread " << count_covered << " samples, total samples " << LandmarkSet.size() << endl;
	cout << "Total boundary parts: " << NumOfBoundaryPart << endl;
	if(count_covered != LandmarkSet.size())
	{
		cout << "Spread part ID wrong" << endl;
		exit(0);
	}
	for (int i = 0; i < nSensors; ++i)
	{
		sensorPool[i].broadcastCovered = 0;
	}
}