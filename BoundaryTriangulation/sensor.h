#ifndef SENSOR_H
#define SENSOR_H

#pragma warning(disable: 4786)

#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <fstream>
#include "math.h"
#include "geometry3D.h"
#include "common.h"

using namespace std;

class Sensor;
class Field;
class LineSegment;
class Triangle;
class Hole;
class Path;
class CircleArea;


typedef map<int, pair<int, int> > DEST2SP; //1st int: dest; 2nd: path length; 3rd: next hop
typedef pair<int, int> SPLEN_NEXTHOP;

typedef set<Sensor*>::const_iterator SCITER;
typedef set<int>::const_iterator IntSetPos;
typedef vector<Sensor*> VSensor;

class Sensor
{
public:
	Sensor();
	~Sensor();

	int index;
	Field *field;
	Point location; //real location
	bool Canceled;
	bool LineSegmentNeighborState;		//used to detect holes
	vector<int> BelongInHoleID;

	set<int>levelsBelongTo;
	bool isBelongToType(int itype);
	int inCommonType(Sensor* s);
	map<int,bool>iType2Landmark;
	map<int,bool>iType2inactive;
	//map<int,int> iType2broadcastCovered;
	map<int,set<int> >iType2VoronoiCellID;

	map<int,int>iType2LandmarkIndex;
	map<int, set<Sensor*> > iTypeLandmarkNeighbors;
	map<int,set<LineSegment*> >iTypeAttachedEdges;
	map<int,set<Triangle* > >iTypeAttachedTriangles;
	map<int, set<Sensor*> > iTypeSickLandmarks;

	void floodingOnItypeBoundaryRouteTo(const int itype, Sensor* dest, vector<Sensor*>& boundPath);
	void floodingRouteOnFullTopology(Sensor* destSensor, vector<Sensor*>& hops);
	void flooding2findAchild();
	void findAChild();
	vector<Sensor*> toChildPath;

	int floodingRouteOnFullTopologyLength(Sensor* destSensor);
	
	Sensor* downLevelLm;
	Sensor* upLevelLm;

	set<Sensor*> neighbors;
	int convex_neighbor_size;
	multimap<double, Sensor*, greater<double> > OnEdgeNeighbors;
	vector<Sensor*> Corresp3rdPoint;				//used to cancel intersect line

	set<Sensor*> KHopNeighbors;
	set<Sensor*> TotalK_1HopNeighbors;
	set<Sensor*> LandmarkNeighbors;

	set<Sensor*> HoleSetNeighbors;

	set<LineSegment*> LandmarkVicilityLines;

	double onEdge_Degree; //on boundary

	void sphericalRoutingOnSameType(Sensor* destSensor,vector<Sensor*>&sameTypePath);
	void sphericalRouting(Sensor* destSensor, vector<Sensor*>& routingPath);
	
	void treeRoutingToDest(Sensor* destSensor, vector<Sensor*>& treePath);
	int Sensor::distanceOnTree(Sensor* destSensor);
	Sensor* treeNexthop(Sensor* destSensor);

	void generateSpanningtree();
	Sensor* parentOnTree;
	bool isChildOnTree;
	set<Sensor*> childrenOnTree;

	void namingSensors();
	vector<Sensor*>name;

	void iTypeLandmarkRouteToLandmark(const int itype,Sensor* destLandmark,vector<Sensor*>& landmarkRoutingPath);
	Sensor* iTypeLandmarkGreedyNextLandmark(const int itype,Sensor* destLandmak);
	long double iTypeSphericalDistance(const int itype, Sensor* destSensor);
	//Added here, which boundary is the sensor nodes supposed to be?
	// -1 not boundary;
	// 0 out boundary;
	// > 0 inner boundary;
	int boundarySensorType;
	int landmarkIndex;

	map<int,DEST2SP>spMapOnItypeBoundary;
	map<int,DEST2SP>spMapOnItypeValidLandmarks;
	map<int,DEST2SP>spMapOnItypeEdgeAttachedLandmarks;

	DEST2SP spMapOnBoundary;
	DEST2SP spMapOnFullTopo;
	DEST2SP spMapOnValidLandmarks;

	DEST2SP spMapOnEdgeAttachedLandmarks;
	int floodOnEdgeAttachedLandmarkNbs(Sensor* dest);
	int floodOnItypeEdgeAttachedLandmarkNbs(const int itype, Sensor* dest);
	void floodOnEdgeAttachedLandmarkNbsRoutes(Sensor* dest, vector<Sensor*>& hopPath);
	void floodOnItypeEdgeAttachedLandmarkNbsRoutes(const int itype, Sensor* dest,vector<Sensor*>& hopPath);

	set<LineSegment*> attachedEdges;
	set<Triangle*> attachedTriangles;
	void find2AccessedSensorInAttachedTriangle(Sensor* &s1, Sensor* &s2);
	LineSegment* findAnAccessedEdgeInAttachedTriangle();

	LineSegment* findAnItypeAccessedEdgeInAttachedTriangle(const int itype);
	void calPosWithCircleIntersecCircle(const Point& c1,const double r1,const Point& c2,const double r2);
	void calPosWithCircleIntersecCircleItype(const int itype,const Point& c1,const double r1,const Point& c2,const double r2);
	set<Sensor*>sickLandmarks;
	int oppositeSideLength;

	
	map<int, long double >iTypeTargetCurvature;
	map<int, long double >iTypeCurrGCurvature;
	map<int, long double >iTypeTheta;
	map<int, Point>iTypeRicciLocation;
	map<int, Point>iTypeSphereLocation;
	map<int, Point>iTypeConformalLocation;
	map<int, long double> iTypeU;
	

	long double targetCurvature;
	long double currGCurvature;
	long double u;
	long double theta;
	bool accessed;

	Point ricciLocation;
	Point sphereLocation;
	// Point unitDiskLocation;
	Point conformalLocation;

	set<Sensor*>edgeAttachedLandmarkNb;
	map<int, set<Sensor*> >itypeEdgeAttachedLandmarkNb;

	void iTypeConforma2Sphere(const int itype);

	void conformal2Sphere();

	void sphere2Unitdisk();
	//In the landmark selection, the important node flood at the itype boundary which is formed by flooding.(second para)
	void floodingOnBoundary(const int itype, const int kHOP, vector<Sensor*>&);

	void floodingOnBoundary(const int itype, const int kHOP);//Flooding khops on the boundary.
	void floodingOnSickLandmarks(Sensor* dest, vector<Sensor*>& holePath);

	void floodingOnItypeSickLandmarks(const int itype, Sensor* dest, vector<Sensor*>& holePath);

	int distanceBetweenLandmarks(const int itype,int khop);// Flooding on the boundary to see if distance < khop.
	void constructVoroniTile(const int itype,int offset );
	// Routing on the boundary, and return the path.
	void floodingOnBoundaryRouteTo(Sensor* dest, vector<Sensor*>& boundPath);


	//temporary variables specific to a flooding
	int broadcastCovered; 
	Sensor *parent;
	Sensor *hole_parent;
	set<Sensor*> children;
	set<Sensor*> candidateParents;

	vector<int> KnownHopDistEndID;//just for preventing surplus computing
	map<int, int> ID_HopDist_List;
	multimap<int, int, less<int> > DistanceToLandmarks; //multimap<Hops, Landmark_ID>

	int level;
	int AvghopToCurrentline;
	int khopConvexNeighborSize;
	multimap<int, int> line_ends_ID;
	double Criticality;
	double Rou_p;


	bool landmark;
	bool inactive;
	bool convex;
	bool concave;
	bool saddle;
	bool asCenter;
	bool asSkeleton;
	bool Be_used_in_triangulation;
	bool DismissedByLastPhase;

	bool within_circle;
	bool put_in_pathcandidateset;

	bool on_boundary_circle;
	bool Out;

	//Ellipsoid
	int boundary_partID;
	vector<Sensor*> SpreadPartID(int part_id);
	bool picked;

	set<int> VoronoiCellID;
	set<int> ConnectedID;

	void clear();
	void init(int index_, Point& location_);
	void updateNeighbors(char* line);
	void updateNeighborsQuasiUBG();

	bool CanFindASensorInLandmarkNeighbors(Sensor *s);
	bool CanFindASensorLandmarkNeighbors(Sensor *t);
	bool DismissThisSensor(vector<Triangle> _Triangle);
	bool DismissThisSensorByHop(LineSegment *currentline, int offset);
	//flooding
	void HoleSensorFlooding(set<Sensor*> hole_sensor_set);
	void LocalFlooding(int khop);
	void LocalFloodingforSkeleton(int khop);
	int LocalFloodingForRefineBoundaryPoints(int khop);
	int LocalFloodingForRefineConvexPoints(int khop);
	int LocalFloodingForRefineConcavePoints(int khop);
	//Sensor* LocalFlooding(int khop);
	


	void LocalFloodingForEstablishCell(int offset);
	void LocalFloodingForEstablishCellforSkeleton(int offset);
	int LocalFloodingforGettingOffsetDistance(Sensor *start, Sensor *end);
	void LocalFloodingForMeasureTheHopsToSensors(int khop);

	vector<Sensor*> PathFlood(Sensor *end);
	vector<Sensor*> PathFloodforDirectPath(Sensor *end);
	vector<Sensor*> PathCriticalityFlooding(Sensor *end, set<Sensor*> area_sensor_set);
	

	int LandmarkFlooding(int khop);
	int FloodingforGetHopsBetweenStartandEnd(Sensor *end);
	int FloodingforGetHopsBetweenStartandEndforSkeleton(Sensor *end);

	bool NotNearEndsOfPath(Path *ph);
	bool SearchInKhopNeighbors(Path *original, Path *pa, int hopnum);

	bool WithinArea(set<Sensor*> area_sensor_set);
	bool PathByBoundarySensors(Sensor *end);

	bool UnknowDistTo(Sensor *end);		//just for computing avg hop distance in circle

	//detect saddle points
	bool FloodingtoFindConcavepointsWithin(int dist, set<Sensor*> &concaveset);

	//get neighborhood
	void getMultiHopNeighborhood(int maxhopCount, set<Sensor*> & nHopNeighbors); //not including the caller himself
	int  getMultiHopNeighborhoodSize(int hopCount);
	void getThe_Kth_HopNeighborhood(int KHopCount);
	void KHOPConvexNeighborSize(int khop);
	bool NotConnectWith(Sensor *end);
	bool BelongOtherIntersectLine(int ID1, int ID2);
	
	inline bool operator == (const Sensor &right) const
	{
		if(index == right.index)
			return true;
		else
			return false;
	}


};

class Path
{
public:
	Path() {path.clear(); start = NULL; end = NULL; pafield = NULL;}
	Path(Sensor *start_, Sensor *end_) {path.clear(); start = start_; end = end_;}
	vector<Sensor*> path;
	Sensor *start;
	Sensor *end;
	Field *pafield;
	void EstablishPath(Sensor *s1, Sensor *s2);
	void EstablishDirectPath(Sensor *s1, Sensor *s2);
	void EstablishPathBySideStep(Sensor *s1, Sensor *s2, set<Sensor*> area_sensor_set);
	void EstablishLowestCriticalityShortestPath(Sensor *_start, Sensor *_end, set<Sensor*> area_sensor_set);
	bool Path_Changed(Sensor *floodsensor, set<Sensor*> area_sensor_set);
	bool PathCrossing(int hopnum);
	bool IntersectWith(LineSegment *dl, int hopnum);
};

class LineSegment
{
public:
	LineSegment(int LineID_, Sensor *start_, Sensor *end_);
	LineSegment();
	~LineSegment();
	// Added 2012-6-28 21:21:38
	vector<Sensor*>boundaryPath;
	void floodStart2EndPathOnBoundary();

	map<int, vector<Sensor*> >itypeBoundaryPath;


	void floodStart2EndPathOnItypeBoundary(const int itype);
	bool valid;

	map<int,long double>iTypeFlowLength;

	long double flowLength;
    long double lineImportance;
	int broadcastedCovered;
	bool atBoundary;

	Sensor* CClockStart;
	Sensor* CClockEnd;

	bool isSharingAVertexWith(LineSegment& otherL);
	Sensor* sharedVertice(LineSegment& otherL);
	set<Triangle*> attachedTriangles;

	map<int, set<Triangle*> >itypeEdgeAttachedTriangles;

	LineSegment(const LineSegment &);

	int LineID;
	Sensor *start;
	Sensor *end;
	Path Lpath;
	Field *line_Field;// in the global field.
	vector<Sensor*> The3rdPoint;

	void init(int LineID_, Sensor *start_, Sensor *end_);
	void ProcessingThe3rdPoint();

//	bool The3rdPointIsValid(Sensor *s, bool canfindgood3rdvertex);
	bool CanFindASensorInThe3rdPoint(Sensor *s);
	Sensor* GetTheGreatest3rdPoint(multimap< double, Sensor*, greater<double> > sharepoints);
	Sensor* GetTheGreatest3rdPointforSkeleton(multimap< double, Sensor*, greater<double> > sharepoints);

	int usedtimes;
	int Processingtimes;
	bool ProcessingLater;
	bool RefinePhaseCanceled;

	double distanceToPoint(Point pt);

	bool containEnd(Point pt); //pt is one of the ends

	bool TheSameLine(LineSegment *l2);
	bool CorrespThe3rdPointOfLineVertex(Sensor *s);
	
	
	inline bool operator == (const LineSegment &ls) 
	{
		return ( (start->index == ls.start->index && end->index == ls.end->index) || (start->index == ls.end->index && end->index == ls.start->index) );
	}
};

class Field
{
public:
	Field();
	Field(Field &field);
	~Field();					//communication range
	
	int MAX_BOUNDARY_TYPE; //
	int ROOT;
	int CORE_RADIUS;
	int banLevel;
	int STARTLEVEL;

	void findAChildForSensors();
	
	vector<Sensor*>removedSurfaceNodes;

	map<int,Sensor* >iTypeRemovedSurfaceNodes;
	map<int,vector<Sensor*> > iType2LandmarkSet;
	map<int,vector<LineSegment> > iTypeEdgesSet;
	map<int,vector<Triangle> > iTypeTriangleSet;
	map<int,vector<LineSegment> >iTypeNoCrossEdgePool;

	map<int, set<int> > iTypeHoleset;// 3hole in the itype triangle set
	void findItype3Hole(const int itype);
	void connectHole2SickLm(const int itype);
	map<int,vector<Triangle> >iTypeRemovedTris;
	// vector<LineSegment*> LandmarkEdges;
	bool isLeadingToIntersectionItype(const int itype, LineSegment& tmpLine);
	bool isEveryItypeEdgeSaturated(const int itype);

	int findEdgeInItypeNoCrossPool(const int itype, Sensor* start_, Sensor* end_);
	void removeItypeTeras(const int itype);





	void initDeployment(const char* topoFile1, const char* topoFile2);
	void calSensorImportance();
	void selectLandMarks(const int kHOP);
	Sensor* ReselectingFloodingSensor(vector<Sensor*>& iTypeBoundarySensors);

	Sensor* findMostImportantNode(const int iType, vector<Sensor*> &);
	void egdesOfLandmarks();
	void landmarkTriangulation();


	void dumpLayersOfTriangle(string file_out_name);
	bool isLeadingToIntersection(LineSegment& tmpLine);
	bool isEveryEdgeSaturated();

	void ricciFlow();
	void dumpRicciLocation(string file_out_name);
	

	void iTypeRicciFlow();
	void dumpItypeRicciLocation(string file_out_name);
	
	void stereographicProjection();
	void ItypeStereographicProjection();

	double tryTreeRoute(int src, int dest, bool dump = false);
	double trySphereRoute(int src,int dest,bool dump = false);
	void dumpShortestandRealPath(Sensor* src,vector<Sensor*> & realPath,vector<Sensor*> & shortestPath,char* filename );

	void dumpSphereLocation(string file_out_name);
	void dumpItypeSphereLocation(string file_out_name);

	vector<LineSegment>noCrossEdgePool;
	void buildNoCrossingEdgePool();
	void dumpTriangleWithID(char * filename);
	void dumpEdges(char * filename);

	void updateTriInfo();

	void testingTriangleOutPut();
	int test3Edges1usedtimeHole();

	int test3ItypeEdges1usedtimeHole(const int itype);

	void itypeAttachTriangleToEachLandmark(const int itype);
	void itypeAttachEdgesToEachLandmark(const int itype);
	void itypeAttachEdgesToTriangles(const int itype);

	int findItypeEdge(const int itype,Sensor* start_, Sensor* end_);
	LineSegment* findItypeEdgePointer(const int itype,Sensor* start_, Sensor* end_);
	int findItypeTriangle(const int itype,Sensor* t1_, Sensor* t2_, Sensor* t3_);
	Triangle* findItypeTrianglePointer(const int itype,Sensor* t1_, Sensor* t2_, Sensor* t3_);

	void sphereAlignment();

	void getUpLowLevelNbLandmarks();




	void attachTriangleToEachLandmark();
	void attachEdgesToEachLandmark();
	void attachEdgesToTriangles();

	void attachEachSensorWithSickLandmarks();
	void itypeAttachattachEachSensorWithSickLandmarks(const int itype);

	void buildEdgeAttachedLandmarkNb();
	void buildItypeEdgeAttachedLandmarkNb(const int itype);
	//void 
	
	// i type landmarks
	

	void removeTeras();
	vector<Triangle> removedTris;
	// remove hole size with 3 by remove a triangle.
	void removeImpossibleHoles();

	int findEdge(Sensor* start_, Sensor* end_);
	LineSegment* findEdgePointer(Sensor* start_, Sensor* end_);
	int findTriangle(Sensor* t1_, Sensor* t2_, Sensor* t3_);
	Triangle* findTrianglePointer(Sensor* t1_, Sensor* t2_, Sensor* t3_);

	int findEdgeInNoCrossPool(Sensor* start_, Sensor* end_);
	// int findTriangleInNoCrossPool(Sensor* t1_, Sensor* t2_, Sensor* t3_);
	// Find the edge in the edge pool
	// bool isLegalEdge(LineSegment&);


	void determineEdgeNodes();
	void refineEdgeNodes();

	void dumpTopology(char * filename);
	void dumpTopologyFull(char * filename);
	void recoverPreprocessedNetwork(char* filename); //read the network from a file
	void recoverPreprocessedTriangle(char *filename);
	void dumpTriangle(char * filename);
	void dumpPath(char * filename, Path path);
	void dumpCircle(char * filename, CircleArea ca);
	void dumpLine(char *filename, Sensor *start, Sensor *end);

	void SelectingLamdmarksRandomly(int khop);
	void SelectingLamdmarksRandomlyForSkeleton(int khop);
	void PickingLandmarkWithSignificance(int khop);
	void detectConcavePoint(int rhop, int dealingpoints);
	void refineConvexpoints(int khop, int within_num);
	void refineConcavepoints(int khop, int within_num);
	void detectSaddlepoints(int detect_dist, int radius);
	void ComputeFieldAvghopDistWithinCircle(int rhop);

	void GenerateVoronoiCell();
	void GenerateVoronoiCellForSkeleton();
	void RefineSurface();
	void RefineSurfaceForSkeleton();

	void EstablishLandmarkVicilitylines();

	void GetRefineSensorAndEdgeSet();

	void GetRefineEdgeSet();
	void DetectHoles();
	void DetectSkeleton(int khop, int offset);

	void ReportTopologyState();
	void DetectEdgeStateforSkeleton();
	//void UpdateRefineEdgeSet();
	void UpdateRefineEdgeSet(vector<LineSegment*> HoleEdgeSet);
	void ClassifyHoleSensors();
	void FillHoles();

	void ResetCancelState();
	void Triangulation();
	void TriangulationForSkeleton();

	bool EndGenerating();
	bool EndUpdating();
	bool LandmarkBelongToRefineSensorSet(Sensor *l);
	bool noCrossingBy(LineSegment *baseline, Sensor* s);
	bool noCrossingByforSkeleton(LineSegment *baseline, Sensor* s);
	bool NotLandmarkNeighbors(int id1, int id2);
	bool PathCrossingInnerNode(int id1, int id2);
	bool WithinTopology(Point Ore);

	set<LineSegment*> TriangleEnclosed(LineSegment *l);
	//Ellipsoid processing
	void PartitionBoundarySamples();
	int NumOfBoundaryPart;

	Sensor* ReselectingFloodingSensor();
	

	int nSensors;
	double comRadius;
	Sensor *sensorPool;
	vector<Sensor*> BoundarySensorSet;
	vector<Sensor*> LandmarkSet;
	multimap<int, int> ConnectingLandmarkIDPairs;

	vector<LineSegment> EdgeSet;
	vector<LineSegment> FinalEdgeSet;
	vector<Triangle> TriangleSet;
	vector<Hole> HoleSet;

	vector<LineSegment*> RefineEdgeSet;
	set<Sensor*> RefineSensorSet;
	set<Sensor*> VoronoiCellSensorSet;
	vector<Sensor*> SkeletonNodeSet;

	double avg1HopNeighborhoodSize;
	double avgNHopNeighborhoodSize;
	double avgOverallNodeDegree;
	double avgOverallHopDistWithinCircle;

	double Avg_distance_err;
	int num_of_outer_boundary_surface;

};


class Triangle
{
public:
	Triangle(int _Tri_ID, Sensor *t1_, Sensor *t2_, Sensor *t3_);
	Triangle();
	~Triangle();
	int Tri_ID;
	Field *Tri_field;
	Sensor *t1;
	Sensor *t2;
	Sensor *t3;
	double Normalx;
	double Normaly;
	double Normalz;

	map<Sensor*,Sensor*>v2vMap;
	map<Sensor*,long double>v2theta;

	map<int,map<Sensor*, Sensor*> >iTypeV2Vmap;
	map<int,map<Sensor*,long double> >iTypeV2Theta;

	map<int, set<Triangle*> >iTypeTriangleNeighbors;
	Triangle* findAItypeNeighborInFrontier(const int itype, vector<Triangle*> &frontier);
	Triangle* findANeighborInFrontier(vector<Triangle*> &frontier);
	LineSegment* findASharingEdgeWith(Triangle* T);

	LineSegment* findASharingItypeEdgeWith(const int itype, Triangle* T);
	Sensor* get3rdSensorByEdge(LineSegment* oppEdge);
	double triImportance;
	bool broadcastCovered;

	void iTypeEstablishTriNeighborhood(const int itype);
	
	bool ContainSensor(Sensor *s);
	void EstablishNeighborhood();
	void GetNormalVector();
	set<Triangle*> TriangleNeighbors;
	inline bool operator == (const Triangle &ls) 
	{
		return ( (t1->index == ls.t1->index && t2->index == ls.t2->index && t3->index == ls.t3->index) 
			|| (t1->index == ls.t1->index && t2->index == ls.t3->index && t3->index == ls.t2->index) 
			|| (t1->index == ls.t2->index && t2->index == ls.t1->index && t3->index == ls.t3->index) 
			|| (t1->index == ls.t2->index && t2->index == ls.t3->index && t3->index == ls.t1->index)
			|| (t1->index == ls.t3->index && t2->index == ls.t1->index && t3->index == ls.t2->index)
			|| (t1->index == ls.t3->index && t2->index == ls.t2->index && t3->index == ls.t1->index) );
	}
	inline bool operator != (const Triangle &ls) 
	{
		return !( (t1->index == ls.t1->index && t2->index == ls.t2->index && t3->index == ls.t3->index) 
			|| (t1->index == ls.t1->index && t2->index == ls.t3->index && t3->index == ls.t2->index) 
			|| (t1->index == ls.t2->index && t2->index == ls.t1->index && t3->index == ls.t3->index) 
			|| (t1->index == ls.t2->index && t2->index == ls.t3->index && t3->index == ls.t1->index)
			|| (t1->index == ls.t3->index && t2->index == ls.t1->index && t3->index == ls.t2->index)
			|| (t1->index == ls.t3->index && t2->index == ls.t2->index && t3->index == ls.t1->index) );
	}
};

class Hole
{
public:
	Hole(int holeID_, set<Sensor*> HoleSensorSet_);
	Hole();
	~Hole();
	int holeID;
	set<Sensor*> HoleSensorSet;
};

class CircleArea
{
public:
	CircleArea();
	CircleArea(Field *_field, Sensor *_center, int _radiu);
	~CircleArea();
	int Radiu;
	double AvgCriticality;
	double AvgHopdist;
	Sensor *Center;
	Field *cirfield;
	int contain_convex_number;
	set<Sensor*> AreaSensorSet;
	multimap<double, Sensor*, less<double> > Circle;
	multimap<double, Sensor*, less<double> > GetCircle();
	void ComputingCircleAvgHopDist();
};




extern Field *globalField;
extern double Low_Edge_Threshold;
extern double High_Edge_Threshold;
extern double Convex_threshold;
extern double Concave_threshold;
extern double Concave_threshold1;
extern double Concave_threshold2;
extern double Saddle_threshold;

extern int KHOP;
extern int Total_Number_of_Triangles;
extern int RemainingEdges;
extern int Number_of_Landmark;

extern int LENGTH;

//debugging states
//extern bool CannotFindASensorInLandmarkNeighbors;
//extern bool CrossingOtherLines;
//extern bool DismissThe3rdPoint;
extern int Contain_convex_threshold;
extern double T_NS;
extern double Gridinterval;
extern double RATIO;
extern FILE *g_fp;

#endif