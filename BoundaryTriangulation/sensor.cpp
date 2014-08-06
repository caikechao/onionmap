#include "sensor.h"
#include "random.h"
#include "assert.h"
#include "common.h"
#include "geometry3D.h"

Sensor::Sensor()
{	
	boundarySensorType = -1;
	// active = false;
	// -1 not boundary;
	// 0 out boundary;
	// > 0 inner boundary;
	landmarkIndex = -1;
	attachedEdges.clear();
	attachedTriangles.clear();
	sickLandmarks.clear();

	spMapOnBoundary.clear();
	spMapOnValidLandmarks.clear();
	oppositeSideLength = 0;
	targetCurvature = 0.0;
	currGCurvature = 0.0;
	u = 0.0;
	theta = 0.0;
	accessed = false;
	//-1 means undefined.

	Canceled = false;
	LineSegmentNeighborState = false;
	Out = false;

	landmark = false;
	inactive = false;
	convex = false;
	concave = false;
	saddle = false;
	within_circle = false;
	put_in_pathcandidateset = false;
	on_boundary_circle = false;
	asCenter = false;
	asSkeleton = false;
	Be_used_in_triangulation = false;

	picked = false;

	index = -1;							
	location = INVALID_POINT;
	neighbors.clear();				
	children.clear();
	LandmarkNeighbors.clear();
	BelongInHoleID.clear();
	parent = NULL;
	hole_parent = NULL;
	broadcastCovered = 0;
	onEdge_Degree = 0;
	Criticality = 1;
	Rou_p = 1;
	AvghopToCurrentline = 0;
	convex_neighbor_size = 0;
	khopConvexNeighborSize = 0;
	boundary_partID = 0;

	OnEdgeNeighbors.clear();
	KHopNeighbors.clear();
	TotalK_1HopNeighbors.clear();
	VoronoiCellID.clear();
	LandmarkVicilityLines.clear();
	KnownHopDistEndID.clear();
	ID_HopDist_List.clear();
	line_ends_ID.clear();
	DistanceToLandmarks.clear();
	ConnectedID.clear();
	//////////////////////////////////////////////////////////////////////////
	//cannotfindsensorsinLandmarkNeighbors = false;
	//crossingOtherLines = false;
	//dismissThe3rdPoint = false;

	parentOnTree = NULL;
	childrenOnTree.clear();
	isChildOnTree = false;
	name.clear();

	downLevelLm = NULL;// next level,level is from root to leaf, the level becomes larger.
	upLevelLm = NULL;
}

Sensor::~Sensor()
{
}
void Sensor::calPosWithCircleIntersecCircleItype(const int itype,const Point& c0,const double r0,const Point& c1,const double r1)
{
	/************************************************************************/
	/* Here ,  the sequence of c0 and c1 is very important...
	very important... it is the key to choose the right solution.!!!!!!!!!
	*/
	/************************************************************************/
	long double a, dx, dy, d, h, rx, ry;
    long double x2, y2;

  /* dx and dy are the vertical and horizontal distances between
   * the circle centers.
   */
  dx = c1.x - c0.x;
  dy = c1.y - c0.y;

  /* Determine the straight-line distance between the centers. */
  //d = sqrt((dy*dy) + (dx*dx));
  d = hypot(dx,dy); // Suggested by Keith Briggs

  /* Check for solvability. */
  if (d > (r0 + r1))
  {
    /* no solution. circles do not intersect. */
	  cout<<"The circle do not intersect!"<<endl;
	  cout<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
	  exit(EXIT_FAILURE);
  }
  if (d < fabs(r0 - r1))
  {
	  cout<<"/* no solution. one circle is contained in the other */"<<endl;
	  cout<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
	  exit(EXIT_FAILURE);
  }

  /* 'point 2' is the point where the line through the circle
   * intersection points crosses the line between the circle
   * centers.  
   */

  /* Determine the distance from point 0 to point 2. */
  a = ((r0*r0) - (r1*r1) + (d*d)) / (2.0 * d) ;

  /* Determine the coordinates of point 2. */
  x2 = c0.x + (dx * a/d);
  y2 = c0.y + (dy * a/d);

  /* Determine the distance from point 2 to either of the
   * intersection points.
   */
  h = sqrt((r0*r0) - (a*a));

  /* Now determine the offsets of the intersection points from
   * point 2.
   */
  rx = -dy * (h/d);
  ry = dx * (h/d);

  /* Determine the absolute intersection points. */
  long double xi = x2 + rx;
  long double xi_prime = x2 - rx;
  long double yi = y2 + ry;
  long double yi_prime = y2 - ry;

  Point i1(xi,yi);
  Point i2(xi_prime,yi_prime);

  /************************************************************************/
  /* there are two points, however, the right solution should satisfy the condition below:

  (the input is v0,v_,v1),then (v_-v0)x(v1-v0)>0
  */
  /************************************************************************/
  Point judgeS = (i1 - c0).crossProduct(c1 - c0);

  if (judgeS.z > 0)
  {
	  iTypeRicciLocation[itype] = i1;
  }
  else
	   iTypeRicciLocation[itype] = i2;

}
void Sensor::calPosWithCircleIntersecCircle(const Point& c0,const double r0,const Point& c1,const double r1)
{
	/************************************************************************/
	/* Here ,  the sequence of c0 and c1 is very important...
	very important... it is the key to choose the right solution.!!!!!!!!!
	*/
	/************************************************************************/
	long double a, dx, dy, d, h, rx, ry;
    long double x2, y2;

  /* dx and dy are the vertical and horizontal distances between
   * the circle centers.
   */
  dx = c1.x - c0.x;
  dy = c1.y - c0.y;

  /* Determine the straight-line distance between the centers. */
  //d = sqrt((dy*dy) + (dx*dx));
  d = hypot(dx,dy); // Suggested by Keith Briggs

  /* Check for solvability. */
  if (d > (r0 + r1))
  {
    /* no solution. circles do not intersect. */
	  cout<<"The circle do not intersect!"<<endl;
	  cout<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
	  exit(EXIT_FAILURE);
  }
  if (d < fabs(r0 - r1))
  {
	  cout<<"/* no solution. one circle is contained in the other */"<<endl;
	  cout<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
	  exit(EXIT_FAILURE);
  }

  /* 'point 2' is the point where the line through the circle
   * intersection points crosses the line between the circle
   * centers.  
   */

  /* Determine the distance from point 0 to point 2. */
  a = ((r0*r0) - (r1*r1) + (d*d)) / (2.0 * d) ;

  /* Determine the coordinates of point 2. */
  x2 = c0.x + (dx * a/d);
  y2 = c0.y + (dy * a/d);

  /* Determine the distance from point 2 to either of the
   * intersection points.
   */
  h = sqrt((r0*r0) - (a*a));

  /* Now determine the offsets of the intersection points from
   * point 2.
   */
  rx = -dy * (h/d);
  ry = dx * (h/d);

  /* Determine the absolute intersection points. */
  long double xi = x2 + rx;
  long double xi_prime = x2 - rx;
  long double yi = y2 + ry;
  long double yi_prime = y2 - ry;

  Point i1(xi,yi);
  Point i2(xi_prime,yi_prime);

  /************************************************************************/
  /* there are two points, however, the right solution should satisfy the condition below:

  (the input is v0,v_,v1),then (v_-v0)x(v1-v0)>0
  */
  /************************************************************************/
  Point judgeS = (i1 - c0).crossProduct(c1 - c0);

  if (judgeS.z > 0)
  {
	  ricciLocation = i1;
  }
  else
	  ricciLocation = i2;

}
LineSegment* Sensor::findAnItypeAccessedEdgeInAttachedTriangle(const int itype)
{
	LineSegment* EdgePos = NULL;
	if (iTypeAttachedTriangles[itype].empty())
	{
		cerr<<"The attached triangles should not be empty!"<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}

	for (auto itTr = iTypeAttachedTriangles[itype].begin(); itTr!= iTypeAttachedTriangles[itype].end();++itTr)
	{
		Triangle* tr = *itTr;
		set<Sensor*>ThreeV;
		ThreeV.insert(tr->t1);
		ThreeV.insert(tr->t2);
		ThreeV.insert(tr->t3);

		ThreeV.erase(this);

		auto countAccessed = count_if(ThreeV.begin(),ThreeV.end(),[](Sensor* s)->bool
		{
			return (s->accessed);
		});

		if (countAccessed == 2)
		{
			Sensor* s1 = *(ThreeV.begin());
			Sensor* s2 = *(++ThreeV.begin());
			EdgePos = this->field->findItypeEdgePointer(itype,s1,s2);
			break;
		}
		else 
			continue;
	}

	return EdgePos;


}
LineSegment* Sensor::findAnAccessedEdgeInAttachedTriangle()
{
	LineSegment* EdgePos = NULL;
	if (attachedTriangles.empty())
	{
		cerr<<"The attached triangles should not be empty!"<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}

	for (auto itTr = attachedTriangles.begin(); itTr!= attachedTriangles.end();++itTr)
	{
		Triangle* tr = *itTr;
		set<Sensor*>ThreeV;
		ThreeV.insert(tr->t1);
		ThreeV.insert(tr->t2);
		ThreeV.insert(tr->t3);

		ThreeV.erase(this);

		auto countAccessed = count_if(ThreeV.begin(),ThreeV.end(),[](Sensor* s)->bool
		{
			return (s->accessed);
		});

		if (countAccessed == 2)
		{
			Sensor* s1 = *(ThreeV.begin());
			Sensor* s2 = *(++ThreeV.begin());
			EdgePos = this->field->findEdgePointer(s1,s2);
			break;
		}
		else 
			continue;
	}

	return EdgePos;

}
void Sensor::find2AccessedSensorInAttachedTriangle(Sensor* &s1, Sensor* &s2)
{
	assert(s1==NULL && s2==NULL);
	
	for (auto itTr = attachedTriangles.begin(); itTr!= attachedTriangles.end();++itTr)
	{
		Triangle* tr = *itTr;
		set<Sensor*>ThreeV;
		ThreeV.insert(tr->t1);
		ThreeV.insert(tr->t2);
		ThreeV.insert(tr->t3);

		ThreeV.erase(this);
		auto countAccessed = count_if(ThreeV.begin(),ThreeV.end(),[](Sensor* s)->bool
		{
			return (s->accessed);
		});

		if (countAccessed == 2)
		{
			s1 = *(ThreeV.begin());
			s2 = *(++ThreeV.begin());
		}
		else 
			continue;
	}
}
void Sensor::clear()
{			
	Canceled = false;
	LineSegmentNeighborState = false;

	landmark = false;
	inactive = false;
	convex = false;
	concave = false;
	saddle = false;
	within_circle = false;
	put_in_pathcandidateset = false;
	on_boundary_circle = false;
	asCenter = false;
	Be_used_in_triangulation = false;

	children.clear();
	candidateParents.clear();
	LandmarkNeighbors.clear();
	BelongInHoleID.clear();
	Corresp3rdPoint.clear();
	HoleSetNeighbors.clear();
	parent = NULL;
	broadcastCovered = 0;

	AvghopToCurrentline = 0;
	convex_neighbor_size = 0;
	khopConvexNeighborSize = 0;

	OnEdgeNeighbors.clear();
	KHopNeighbors.clear();
	TotalK_1HopNeighbors.clear();
	VoronoiCellID.clear();
	LandmarkVicilityLines.clear();
	KnownHopDistEndID.clear();
	ID_HopDist_List.clear();
	line_ends_ID.clear();
	DistanceToLandmarks.clear();
	ConnectedID.clear();
}

#if 0
void Sensor::updateNeighbors()
{
	neighbors.clear();

	for(int i = 0; i < field->nSensors; i++) {
		Sensor & s = field->sensorPool[i];			//s is equivalent to field->sensorPool[i]
		if (i != index && location.distance(s.location) < COMM_RANGE ){
			neighbors.insert(&s);
		}
	}
	return;
}
#endif

void Sensor::updateNeighbors(char* line)
{
	neighbors.clear();
	int i = 0;
	int Index;
	int count = 0;
	char val[1024];
	while(line[i] != '\n')
	{
		int j = 0;
		while(line[i] != ' ' && line[i] != '\n')
		{
			val[j] = line[i];
			++i;
			++j;
		}
		Index = atoi(val);
		memset(val, 0, sizeof(val));
		++count;
		if(count == 1 && Index != index)
		{
			fprintf_s(stderr, "Information wrong!\n");
			exit(0);
		}
		if(count != 1)
		{
			neighbors.insert(&field->sensorPool[Index]);
		}
		++i;
	}
			
	return;
}

void Sensor::updateNeighborsQuasiUBG()
{
	
	int i;
	neighbors.clear();

	UniformRandVar rand;
	double prob = (2 - QUASI_UBG_ALPHA) / 3.1;

	for(i = 0; i < field->nSensors; i++) {
		Sensor & s = field->sensorPool[i];
		if (i == index) 
			continue;

		double d = location.distance(s.location);

		if (d < (1 - QUASI_UBG_ALPHA) * COMM_RANGE ){
			neighbors.insert(&s);
		}
		else if (d < (1 + QUASI_UBG_ALPHA) * COMM_RANGE) {
			if (rand.value() < prob) {
				neighbors.insert(&s);
			}
		}
	}
}
bool Sensor::isBelongToType(int itype)
{
	bool belongTo = false;
	if (levelsBelongTo.empty())
	{
		cerr<<"The levels belong to is empty!"<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	auto it = find(levelsBelongTo.begin(),levelsBelongTo.end(),itype);

	if (it == levelsBelongTo.end())
		belongTo = false;
	else
		belongTo = true;

	return belongTo;

}
void Sensor::init(int index_, Point& location_)
{
	index = index_;
	location = location_;

	levelsBelongTo.swap(location_.levelTypes);
	boundarySensorType = location_.boundaryType;
}

void Sensor::getMultiHopNeighborhood(int nHops,set<Sensor*> & nHopNeighbors)
{
	assert(nHops >= 1);
	nHopNeighbors.clear();
	// Get the n hop neighbor hood.
	int hopCounter = 0;
	set<Sensor*>touchedNodes;
	if (nHops == 1)
	{
		for (auto neighborIt = neighbors.begin(); neighborIt != neighbors.end(); ++ neighborIt)
		{
			nHopNeighbors.insert(*neighborIt);
		}
	}
	else 
	{
		set<Sensor*> *frontier = new set<Sensor*>;
		this->broadcastCovered = 1;
		frontier->insert(this);
		touchedNodes.insert(this);

		while(hopCounter < nHops && !frontier->empty())
		{
			set<Sensor*> *newFrontier = new set<Sensor*>;
			for(set<Sensor*>::iterator iter = frontier->begin(); iter != frontier->end(); ++iter)
			{
				Sensor *s = (*iter);
				for(auto snpos = s->neighbors.begin(); snpos != s->neighbors.end(); ++snpos)
				{
					Sensor *n = *snpos;
					if(n->broadcastCovered == 0)
					{
						nHopNeighbors.insert(n);
						newFrontier->insert(n);
						n->broadcastCovered = 1;
						touchedNodes.insert(n);

					}
				}
			}
			delete frontier;
			frontier = newFrontier;
			hopCounter++;
		}
		delete frontier;
	}

	// Reinitialize the broadcasting.
	for (auto touchedIt = touchedNodes.begin(); touchedIt != touchedNodes.end(); ++ touchedIt)
	{
		Sensor* s = *touchedIt;
		s->broadcastCovered = 0;
	}
}


int Sensor::getMultiHopNeighborhoodSize(int hopCount)
{
	assert(hopCount >= 1);
	int neighborsSize = 0;
	// nHopNeighbors.clear();
	// Get the n hop neighbor hood.
	int hopCounter = 0;
	set<Sensor*>touchedNodes;
	if (hopCount == 1)
	{
		neighborsSize = neighbors.size();
	}
	else 
	{
		set<Sensor*> *frontier = new set<Sensor*>;
		this->broadcastCovered = 1;
		frontier->insert(this);
		touchedNodes.insert(this);

		while(hopCounter < hopCount && !frontier->empty())
		{
			set<Sensor*> *newFrontier = new set<Sensor*>;
			for(set<Sensor*>::iterator iter = frontier->begin(); iter != frontier->end(); ++iter)
			{
				Sensor *s = (*iter);
				for(auto snpos = s->neighbors.begin(); snpos != s->neighbors.end(); ++snpos)
				{
					Sensor *n = *snpos;
					if(n->broadcastCovered == 0)
					{
						++ neighborsSize;
						newFrontier->insert(n);
						n->broadcastCovered = 1;
						touchedNodes.insert(n);

					}
				}
			}
			delete frontier;
			frontier = newFrontier;
			hopCounter++;
		}
		delete frontier;
	}

	// Reinitialize the broadcasting.
	for (auto touchedIt = touchedNodes.begin(); touchedIt != touchedNodes.end(); ++ touchedIt)
	{
		Sensor* s = *touchedIt;
		s->broadcastCovered = 0;
	}
	
	return neighborsSize;

}

void Sensor::getThe_Kth_HopNeighborhood(int KHopCount)
{
	SCITER neighborpos;
	int hopCount = 1;
//	while(KHopCount)
//	{
	if (KHopCount == 1) 
	{
		KHopNeighbors.insert(neighbors.begin(), neighbors.end());
		return;
	}
	set<Sensor*> *frontier = new set<Sensor*>;
	frontier->insert(this);
	while (hopCount <= KHopCount && !frontier->empty()) 
	{
		set<Sensor*> *newFrontier = new set<Sensor*>;
		for(SCITER iter = frontier->begin(); iter != frontier->end(); ++iter) 
		{
			Sensor *s = (*iter);
			for(neighborpos = s->neighbors.begin(); neighborpos != s->neighbors.end(); neighborpos++) 
			{
				Sensor *n = *neighborpos;	
				if (n != this && n->broadcastCovered == 0) 
				{
					if(hopCount == KHopCount)
					{
						KHopNeighbors.insert(n);
					}
					else
					{
						TotalK_1HopNeighbors.insert(n);
						n->broadcastCovered = 1;
					}
					newFrontier->insert(n);
				}					
			}
		}
		delete frontier;		
		frontier = newFrontier;			
		hopCount++;
	}

	delete frontier;
//		--KHopCount;
//	}
	for(int i = 0; i < field->nSensors; i++) {
		field->sensorPool[i].broadcastCovered = 0;
	}
}

void Sensor::KHOPConvexNeighborSize(int khop)
{
	int depth = 0;
	set<Sensor*> *frontier = new set<Sensor*>;
	broadcastCovered = 1;
	frontier->insert( this );
	while(depth != khop && !frontier->empty())
	{
		set<Sensor*> *newFrontier = new set<Sensor*>;
		for(set<Sensor*>::const_iterator iFron = frontier->begin(); iFron != frontier->end(); ++iFron) 
		{
			Sensor *f = *iFron;

			for(SCITER NeighborPos = f->neighbors.begin(); NeighborPos != f->neighbors.end(); ++NeighborPos)
			{
				Sensor *s = *NeighborPos;
				if(s->broadcastCovered == 0 && s->onEdge_Degree != 0)
				{
					newFrontier->insert( s );
				}
			}
		}
		if(!newFrontier->empty())
		{
			++depth;
		}

		//dealing newFrontier
		set<Sensor*>::iterator newFrontierpos;
		for(newFrontierpos = newFrontier->begin(); newFrontierpos != newFrontier->end(); ++newFrontierpos)
		{
			Sensor *s = *newFrontierpos;
			s->broadcastCovered = 1;
			if(s->convex)
			{
				this->khopConvexNeighborSize++;
			}

		}
		delete frontier;
		frontier = newFrontier;
	}
	delete frontier;
	for(int i = 0; i < field->nSensors; ++i)
	{
		field->sensorPool[i].broadcastCovered = 0;
	}
}

bool Sensor::CanFindASensorInLandmarkNeighbors(Sensor *s)
{
	SCITER LMneighborpos;
	for(LMneighborpos = LandmarkNeighbors.begin(); LMneighborpos != LandmarkNeighbors.end(); ++LMneighborpos)
	{
		Sensor *sen = *LMneighborpos;
		if(sen->index == s->index)
		{
			printf("Can find end's landmark neighbor\n");
			//CannotFindASensorInLandmarkNeighbors = false;
			return true;
		}
	}
	return false;
}

bool Sensor::CanFindASensorLandmarkNeighbors(Sensor *t)
{
	SCITER landmarkneighborpos;
	for(landmarkneighborpos = LandmarkNeighbors.begin(); landmarkneighborpos != LandmarkNeighbors.end(); ++landmarkneighborpos)
	{
		Sensor *sen = *landmarkneighborpos;
		if(sen->index == t->index)
			return true;
	}
	return false;
}

bool Sensor::DismissThisSensor(vector<Triangle> _Triangle)
{
	if(_Triangle.empty() || this->convex)
		return false;
	else
	{
		for(unsigned int i = 0; i < _Triangle.size(); ++i)
		{
			if(_Triangle[i].ContainSensor(this))
			{
				printf("DismissThisSensor\n");
				return true;
			}
		}
	}
	return false;
}
#if 0
bool Sensor::DismissThisSensorByHop(LineSegment *currentline, int offset)
{
	int avghop = 0;

	int _3rdhop = 0;
	bool _3rdhoptag = false;

	int starthop = 0;
	bool starthoptag = false;

	int endhop = 0;
	bool endhoptag = false;

	int depth = 0;
	if(currentline->The3rdPoint.empty())
	{
		DismissThe3rdPoint = false;
		return false;
	}
	else if(currentline->The3rdPoint.size() != 1)
	{
		printf("currentline wrong\n");
		exit(0);
	}
	multimap<double, Sensor*, greater<double> > *frontier = new multimap<double, Sensor*, greater<double> >;
	frontier->insert( make_pair(this->onEdge_Degree, this) );
	while(!frontier->empty())
	{
		multimap<double, Sensor*, greater<double> > *newFrontier = new multimap<double, Sensor*, greater<double> >;
		for(multimap<double, Sensor*, greater<double> >::const_iterator iFron = frontier->begin(); iFron != frontier->end(); ++iFron) 
		{
			Sensor *f = iFron->second;

			for(SCITER NeighborPos = f->neighbors.begin(); NeighborPos != f->neighbors.end(); ++NeighborPos)
			{
				Sensor *s = *NeighborPos;
				if(s->broadcastCovered == 0 && s->onEdge_Degree != 0)
				{
					newFrontier->insert( make_pair(s->onEdge_Degree, s) );
					s->broadcastCovered = 1;
				}
			}
		}
		if(!newFrontier->empty())
		{
			++depth;
		}
		//dealing newFrontier
		multimap<double, Sensor*, greater<double> >::iterator newFrontierpos;
		for(newFrontierpos = newFrontier->begin(); newFrontierpos != newFrontier->end(); ++newFrontierpos)
		{
			Sensor *s = newFrontierpos->second;
			//s->inactive = true;
			s->broadcastCovered = 1;
			if(s->index == currentline->start->index)
			{
				starthop = depth;
				starthoptag = true;
			}
			if(s->index == currentline->end->index)
			{
				endhop = depth;
				endhoptag = true;
			}
			if(s->index == currentline->The3rdPoint[0]->index)
			{
				_3rdhop = depth;
				_3rdhoptag = true;
			}
		}
		if(_3rdhoptag && starthoptag && endhoptag)
		{
			avghop = (starthop + endhop)/2;
			int OFT = abs(avghop - _3rdhop);
			if(_3rdhop < avghop)
			{
				printf("Dismiss the 3rd point\n");
				DismissThe3rdPoint = true;
				return true;
			}
			else if(OFT <= offset)
			{
				printf("Dismiss the 3rd point\n");
				DismissThe3rdPoint = true;
				return true;
			}
			else
			{
				DismissThe3rdPoint = false;
				return false;
			}
		}
		frontier = newFrontier;
	}
	delete frontier;
	for(unsigned int i = 0; i < field->EdgeSensorSet.size(); ++i)
	{
		field->EdgeSensorSet[i]->broadcastCovered = 0;
	}
	DismissThe3rdPoint = false;
	return false;
}
#endif

bool Sensor::DismissThisSensorByHop(LineSegment *currentline, int offset)
{
	int avghop = 0;

	int _3rdhop = 0;
	bool _3rdhoptag = false;

	int starthop = 0;
	bool starthoptag = false;

	int endhop = 0;
	bool endhoptag = false;

	int depth = 0;
	if(currentline->The3rdPoint.empty())
	{
		//DismissThe3rdPoint = false;
		return false;
	}
	else if(currentline->The3rdPoint.size() != 1)
	{
		printf("currentline wrong\n");
		exit(0);
	}
	set<Sensor*> *frontier = new set<Sensor*>;
	
	broadcastCovered = 1;
	frontier->insert( this );
	while(!frontier->empty())
	{
		set<Sensor*> *newFrontier = new set<Sensor*>;
		for(set<Sensor*>::const_iterator iFron = frontier->begin(); iFron != frontier->end(); ++iFron) 
		{
			Sensor *f = *iFron;

			for(SCITER NeighborPos = f->neighbors.begin(); NeighborPos != f->neighbors.end(); ++NeighborPos)
			{
				Sensor *s = *NeighborPos;
				if(s->broadcastCovered == 0)
				{
					newFrontier->insert( s );
				}
			}
		}
		if(!newFrontier->empty())
		{
			++depth;
		}
		//dealing newFrontier
		set<Sensor*>::iterator newFrontierpos;
		for(newFrontierpos = newFrontier->begin(); newFrontierpos != newFrontier->end(); ++newFrontierpos)
		{
			Sensor *s = *newFrontierpos;
			s->broadcastCovered = 1;
			if(s->index == currentline->start->index)
			{
				starthop = depth;
				starthoptag = true;
			}
			if(s->index == currentline->end->index)
			{
				endhop = depth;
				endhoptag = true;
			}
			if(s->index == currentline->The3rdPoint[0]->index)
			{
				_3rdhop = depth;
				_3rdhoptag = true;
			}	
		}
		if(_3rdhoptag && starthoptag && endhoptag)
		{
			avghop = (starthop + endhop)/2;
			this->AvghopToCurrentline = avghop;
			int OFT = abs(avghop - _3rdhop);
			if(_3rdhop < avghop)
			{
				printf("Dismiss the 3rd point\n");
				//DismissThe3rdPoint = true;
				delete newFrontier;
				delete frontier;
				return true;
			}
			/*
			else if(OFT <= offset)
						{
							printf("Dismiss the 3rd point\n");
							DismissThe3rdPoint = true;
							return true;
						}*/
			
			else
			{
				//DismissThe3rdPoint = false;
				delete newFrontier;
				delete frontier;
				return false;
			}
		}
		delete frontier;
		frontier = newFrontier;
	}
	delete frontier;
	for(int i = 0; i < field->nSensors; ++i)
	{
		field->sensorPool[i].broadcastCovered = 0;
	}
	//DismissThe3rdPoint = false;
	return false;
}

//for selecting landmarks version1
#if 0
Sensor* Sensor::LocalFlooding(int khop)
{
	int depth = 0;
	multimap<double, Sensor*, greater<double> > *frontier = new multimap<double, Sensor*, greater<double> >;
	broadcastCovered = 1;
	frontier->insert( make_pair(this->onEdge_Degree, this) );
	while(depth != khop)
	{
		multimap<double, Sensor*, greater<double> > *newFrontier = new multimap<double, Sensor*, greater<double> >;
		for(multimap<double, Sensor*, greater<double> >::const_iterator iFron = frontier->begin(); iFron != frontier->end(); ++iFron) 
		{
			Sensor *f = iFron->second;
			
			for(SCITER NeighborPos = f->neighbors.begin(); NeighborPos != f->neighbors.end(); ++NeighborPos)
			{
				Sensor *s = *NeighborPos;
				if(s->broadcastCovered == 0 && s->onEdge_Degree != 0)
				{
					newFrontier->insert( make_pair(s->onEdge_Degree, s) );
					s->broadcastCovered = 1;
				}
			}
		}
		if(!newFrontier->empty())
		{
			++depth;
		}
		else
			return NULL;
		//dealing newFrontier
		multimap<double, Sensor*, greater<double> >::iterator newFrontierpos;
		for(newFrontierpos = newFrontier->begin(); newFrontierpos != newFrontier->end(); ++newFrontierpos)
		{
			Sensor *s = newFrontierpos->second;
			if(depth == khop)
			{
				s->landmark = true;
				s->broadcastCovered = 1;
				++Number_of_Landmark;
				delete frontier;
				return s;
			}
			else
			{
				s->inactive = true;
				s->broadcastCovered = 1;
			}
		}
		frontier = newFrontier;
	}
	delete frontier;
	return NULL;
}
#endif
void Sensor::floodingOnItypeSickLandmarks(const int itype, Sensor* destSensor, vector<Sensor*>& holePath)
{
	// destSensor should not be null.
	assert(destSensor);
	// hole path should be empty
	assert(holePath.empty());
	// both should be landmarks.
	assert(iType2Landmark[itype] && destSensor->iType2Landmark[itype]);
	
	assert(this->isBelongToType(itype) && destSensor->isBelongToType(itype));

	// int thisBoundaryType = boundarySensorType;
	int i;
	int hopCount = 0;
	set<Sensor*> touchedNodes;
	bool deadLoop = false;
	//trivial case
	if (destSensor == this)
		return;

	vector<Sensor*> *frontier = new vector<Sensor*>;
	frontier->push_back(this);
	touchedNodes.insert(this);
	broadcastCovered = 1;

	while (destSensor->broadcastCovered != 1)
	{
		vector<Sensor*> *newFrontier = new vector<Sensor*>;

		for(i = 0; i < (int)(frontier->size()); i++)
		{
			Sensor *s = (*frontier)[i];

			if (s->broadcastCovered == 1)
			{
				for (auto it_sNb = s->iTypeSickLandmarks[itype].begin(); it_sNb != s->iTypeSickLandmarks[itype].end(); ++ it_sNb)
				{
					Sensor *n = *it_sNb;
					// this condition limited the flooding on the boundary.
					if (n->broadcastCovered == 0 && n->isBelongToType(itype))
					{
						n->spMapOnItypeValidLandmarks[itype].insert(DEST2SP::value_type(index, SPLEN_NEXTHOP(hopCount + 1, s->index)));
						n->broadcastCovered = 1;
						touchedNodes.insert(n);
						newFrontier->push_back(n);
					}
				}
			}
		}

		delete frontier;
		frontier = newFrontier;

		hopCount++;
		if (hopCount > 150)
		{
			deadLoop = true;

			for (auto it = touchedNodes.begin(); it != touchedNodes.end(); ++ it)
			{
				(*it)->broadcastCovered = 0;
			}

			// clear the shortest path map, it can be used only once here.
			for(int i =0; i<field->nSensors; ++i)
			{
				Sensor* s = &(field->sensorPool[i]);
				s->spMapOnItypeValidLandmarks[itype].clear();
			}


			return;
		}
	}

	Sensor* curNode = destSensor;

	do
	{
		holePath.insert(holePath.begin(), curNode);
		Sensor *curParent = &globalField->sensorPool[curNode->spMapOnItypeValidLandmarks[itype][index].second];
		//also remember the reverse route
		curParent->spMapOnItypeValidLandmarks[itype].insert(DEST2SP::value_type(destSensor->index, SPLEN_NEXTHOP(holePath.size(), curNode->index)));
		curNode = curParent;
	}
	while (curNode != this);

	for (auto it = touchedNodes.begin(); it != touchedNodes.end(); ++ it)
	{
		(*it)->broadcastCovered = 0;
	}

	// clear the shortest path map, it can be used only once here.
	for(int i =0; i<field->nSensors; ++i)
	{
		Sensor* s = &(field->sensorPool[i]);
		s->spMapOnItypeValidLandmarks[itype].clear();
	}

}
void Sensor::floodingOnSickLandmarks(Sensor* destSensor, vector<Sensor*>& holePath)
{
	// destSensor should not be null.
	assert(destSensor);
	// hole path should be empty
	assert(holePath.empty());
	// both should be landmarks.
	assert(this->landmark && destSensor->landmark);
	assert(boundarySensorType == destSensor->boundarySensorType);

	int thisBoundaryType = boundarySensorType;
	int i;
	int hopCount = 0;
	set<Sensor*> touchedNodes;

	//trivial case
	if (destSensor == this)
		return;

	vector<Sensor*> *frontier = new vector<Sensor*>;
	frontier->push_back(this);
	touchedNodes.insert(this);
	broadcastCovered = 1;

	while (destSensor->broadcastCovered != 1)
	{
		vector<Sensor*> *newFrontier = new vector<Sensor*>;

		for(i = 0; i < (int)(frontier->size()); i++)
		{
			Sensor *s = (*frontier)[i];

			if (s->broadcastCovered == 1)
			{
				for (auto it_sNb = s->sickLandmarks.begin(); it_sNb != s->sickLandmarks.end(); ++ it_sNb)
				{
					Sensor *n = *it_sNb;
					// this condition limited the flooding on the boundary.
					if (n->broadcastCovered == 0 && thisBoundaryType == n->boundarySensorType)
					{
						n->spMapOnValidLandmarks.insert(DEST2SP::value_type(index, SPLEN_NEXTHOP(hopCount + 1, s->index)));
						n->broadcastCovered = 1;
						touchedNodes.insert(n);
						newFrontier->push_back(n);
					}
				}
			}
		}

		delete frontier;
		frontier = newFrontier;

		hopCount++;
	}

	Sensor* curNode = destSensor;

	do
	{
		holePath.insert(holePath.begin(), curNode);
		Sensor *curParent = &globalField->sensorPool[curNode->spMapOnValidLandmarks[index].second];
		//also remember the reverse route
		curParent->spMapOnValidLandmarks.insert(DEST2SP::value_type(destSensor->index, SPLEN_NEXTHOP(holePath.size(), curNode->index)));
		curNode = curParent;
	}
	while (curNode != this);

	for (auto it = touchedNodes.begin(); it != touchedNodes.end(); ++ it)
	{
		(*it)->broadcastCovered = 0;
	}

	// clear the shortest path map, it can be used only once here.
	for(int i =0; i<field->nSensors; ++i)
	{
		Sensor* s = &(field->sensorPool[i]);
		s->spMapOnValidLandmarks.clear();
	}


}
int Sensor::inCommonType(Sensor* otherS)
{
	// find the inner most common itype surface.
	// -1 means the two sensors are not in the same type.

	int max_set_len = max(levelsBelongTo.size(),otherS->levelsBelongTo.size());
	vector<int>commonTypeSet(max_set_len);
	vector<int>::iterator lastcommonNext;

	lastcommonNext = set_intersection(levelsBelongTo.begin(),levelsBelongTo.end(),
									  otherS->levelsBelongTo.begin(),otherS->levelsBelongTo.end(),commonTypeSet.begin());

	int intersectSize = int(lastcommonNext - commonTypeSet.begin());

	UniformRandVar randv(0, intersectSize - 1);
	int pos = (int)randv.value();

	if (intersectSize > 0)
	{
		return commonTypeSet[pos];
	}
	else
		return -1;

}
Sensor* Sensor::iTypeLandmarkGreedyNextLandmark(const int itype,Sensor* destLandmark)
{
	if (! this->iType2Landmark[itype])
	{
		cerr<<"This function can only be called by a itype landmark."<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}

	long double minDis = iTypeSphereLocation[itype].distance(destLandmark->iTypeSphereLocation[itype]);
	Sensor* nextLm = NULL;

	if (minDis == 0.0)// no two landmarks are in the same position...
	{
		return this;
	}
	for (auto it_nbLM  = itypeEdgeAttachedLandmarkNb[itype].begin(); 
		      it_nbLM != itypeEdgeAttachedLandmarkNb[itype].end(); ++ it_nbLM)
	{
		Sensor* s = *it_nbLM;
		long double d = s->iTypeSphereLocation[itype].distance(destLandmark->iTypeSphereLocation[itype]);

		if (d < minDis)
		{
			minDis = d;
			nextLm = s;
		}

	}

	return nextLm;

}
void Sensor::iTypeLandmarkRouteToLandmark(const int itype,Sensor* destLandmark,vector<Sensor*>& landmarkRoutingPath)
{
	if (! this->iType2Landmark[itype])
	{
		cerr<<"This function can only be called by a itype landmark."<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	// greedy route to the dest landmark.
	long double totalDistance = iTypeSphericalDistance(itype,destLandmark);
	Sensor* currLm = this;
	while (currLm != destLandmark)// must success in the greedy routing...
	{
		Sensor* nextLm = currLm->iTypeLandmarkGreedyNextLandmark(itype,destLandmark);
		vector<Sensor*>middle_path;
		currLm->floodingOnItypeBoundaryRouteTo(itype,nextLm,middle_path);
		landmarkRoutingPath.insert(landmarkRoutingPath.end(),middle_path.begin(),middle_path.end());
		currLm = nextLm;
	}


}
long double Sensor::iTypeSphericalDistance(const int itype, Sensor* destSensor)
{
	if (! this->iType2Landmark[itype] && ! destSensor->iType2Landmark[itype])
	{
		cerr<<"iTypeSphericalDistance function can only be called by a itype landmark."<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	// calculate the spherical distance.
	long double thirdEdgeLen = iTypeSphereLocation[itype].distance(destSensor->iTypeSphereLocation[itype]);
	// as the radius is 1 for a unit sphere, we can calculate the central angle.
	long double centralAngle = acos( (1+1 - pow(thirdEdgeLen,(long double)2.0 ) ) / ( 2*1*1 ));

	// the arc length is ...
	return (1.0 * centralAngle);


}
int Sensor::floodingRouteOnFullTopologyLength(Sensor* destSensor)
{
	int i;
	int hopCount = 0;
	set<int> touchedNodes;
	assert(destSensor);
	vector<Sensor*>hops;
	hops.clear();
	if (destSensor == NULL)
	{
		cerr<<"Error, the dest sensor should not be null."<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}

	hops.clear();

	//trivial case
	if (destSensor == this)
		return 0;


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
				for(auto it_nb = s->neighbors.begin();it_nb != s->neighbors.end(); ++ it_nb)
				{
					Sensor *n = *it_nb;
					if (n->broadcastCovered == 0)
					{
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


	for(set<int>::const_iterator iter = touchedNodes.begin(); iter != touchedNodes.end(); ++iter)
	{
		// the first line is to save memory.
		globalField->sensorPool[*iter].spMapOnFullTopo.clear();
		globalField->sensorPool[*iter].broadcastCovered = 0;
		globalField->sensorPool[*iter].parent = NULL;
	}

	return hopCount;

}
void Sensor::findAChild()
{
	// for a sensor level larger than 1
	int thisLowestLevel = *(levelsBelongTo.begin());
	if (thisLowestLevel > field->MAX_BOUNDARY_TYPE || thisLowestLevel < 0)
	{
		cerr<<"Error, this sensor is at an impossible level."<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	if (thisLowestLevel == field->MAX_BOUNDARY_TYPE)
	{
		// no child.
		toChildPath.clear();
		return;
	}
	if (thisLowestLevel < field->MAX_BOUNDARY_TYPE && thisLowestLevel >= 1)
	{
		// should find a child.
		if (!childrenOnTree.empty())
		{
			toChildPath.push_back(*(childrenOnTree.begin()));
		}
		else
		{
			// finding the path to a larger level child.
			flooding2findAchild();
			if (toChildPath.empty())
			{
				cerr<<"Error, can not find the child of the sensor."<<endl;
				cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
				exit(EXIT_FAILURE);
			}


		}
	}
}
void Sensor::flooding2findAchild()
{
	int i;
	int thisLevel = *(levelsBelongTo.begin());
	int hopCount = 0;
	set<int> touchedNodes;


	toChildPath.clear();


	vector<Sensor*> *frontier = new vector<Sensor*>;
	frontier->push_back(this);
	touchedNodes.insert(index);
	broadcastCovered = 1;
	bool outofLoop = false;
	Sensor* aChild = NULL;
	int nLevel = 0;

	while (true)
	{

		vector<Sensor*> *newFrontier = new vector<Sensor*>;

		for(i = 0; i < (int)(frontier->size()); i++)
		{
			Sensor *s = (*frontier)[i];

			if (s->broadcastCovered == 1)
			{
				for(auto it_nb = s->neighbors.begin();it_nb != s->neighbors.end(); ++ it_nb)
				{
					Sensor *n = *it_nb;


					if (n->broadcastCovered == 0)
					{

						n->spMapOnFullTopo.insert(DEST2SP::value_type(index, SPLEN_NEXTHOP(hopCount + 1, s->index)));
						n->broadcastCovered = 1;
						touchedNodes.insert(n->index);
						newFrontier->push_back(n);
						nLevel = *(n->levelsBelongTo.begin());
						if (nLevel == thisLevel +1 )
						{

							aChild = n;
							outofLoop = true;
							break;

						}

					}
				}
			}

			if (outofLoop)
			{
				break;
			}
		}

		delete frontier;
		frontier = newFrontier;

		hopCount++;

		if (outofLoop)
		{
			break;
		}
	}

	delete frontier;

	Sensor* curNode = aChild;

	do
	{
		toChildPath.insert(toChildPath.begin(), curNode);
		Sensor *curParent = &globalField->sensorPool[curNode->spMapOnFullTopo[index].second];
		//also remember the reverse route
		// curParent->spMapOnFullTopo.insert(DEST2SP::value_type(aChild->index, SPLEN_NEXTHOP(toChildPath.size(), curNode->index)));
		curNode = curParent;
	}
	while (curNode != this);

	for(set<int>::const_iterator iter = touchedNodes.begin(); iter != touchedNodes.end(); ++iter)
	{
		// the first line is to save memory.
		globalField->sensorPool[*iter].spMapOnFullTopo.clear();
		globalField->sensorPool[*iter].broadcastCovered = 0;
		globalField->sensorPool[*iter].parent = NULL;
	}
	
}
void Sensor::floodingRouteOnFullTopology(Sensor* destSensor, vector<Sensor*>& hops)
{
	int i;
	int hopCount = 0;
	set<int> touchedNodes;
	assert(destSensor);
	if (destSensor == NULL)
	{
		cerr<<"Error, the dest sensor should not be null."<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}

	hops.clear();

	//trivial case
	if (destSensor == this)
		return;

	//the path is in the cache
	// to save the memory.
// 	if (spMapOnFullTopo.find(destSensor->index) != spMapOnFullTopo.end())
// 	{
// 		Sensor* curNode = this;
// 
// 		while (curNode != destSensor)
// 		{
// 			Sensor *curParent = &globalField->sensorPool[curNode->spMapOnFullTopo[destSensor->index].second];
// 			hops.push_back(curParent);
// 			curNode = curParent;
// 		};
// 
// 		return;
// 	}

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
				for(auto it_nb = s->neighbors.begin();it_nb != s->neighbors.end(); ++ it_nb)
				{
					Sensor *n = *it_nb;
					

					if (n->broadcastCovered == 0)
					{

						n->spMapOnFullTopo.insert(DEST2SP::value_type(index, SPLEN_NEXTHOP(hopCount + 1, s->index)));
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
		hops.insert(hops.begin(), curNode);
		Sensor *curParent = &globalField->sensorPool[curNode->spMapOnFullTopo[index].second];
		//also remember the reverse route
		// curParent->spMapOnFullTopo.insert(DEST2SP::value_type(destSensor->index, SPLEN_NEXTHOP(hops.size(), curNode->index)));
		curNode = curParent;
	}
	while (curNode != this);

	for(set<int>::const_iterator iter = touchedNodes.begin(); iter != touchedNodes.end(); ++iter)
	{
		// the first line is to save memory.
		globalField->sensorPool[*iter].spMapOnFullTopo.clear();
		globalField->sensorPool[*iter].broadcastCovered = 0;
		globalField->sensorPool[*iter].parent = NULL;
	}
}
void Sensor::generateSpanningtree()
{
	//Find the spanning tree of the sensor graph.
	//The root is the current sensor, should be called by the root sensor.
	// assert(this->isRoot == true);
	for (int i = 0; i < field->nSensors; ++i)
	{
		Sensor* s = field->sensorPool + i;
		s->broadcastCovered = 0;
	}
	Sensor* root =  this;
	set<Sensor*> *frontier = new set<Sensor*>;
	frontier->insert(this);
	broadcastCovered = 1;

	set<Sensor*> touchedSensors;
	parentOnTree = NULL;
	isChildOnTree = false;
	touchedSensors.insert(this);

	while (!frontier->empty())
	{	
		set<Sensor*> *newFrontier = new set<Sensor*>;	
		for (auto it =frontier->begin(); it!=frontier->end(); ++ it)
		{
			Sensor *s = *it;
			for(auto it_nb = s->neighbors.begin();it_nb != s->neighbors.end(); ++ it_nb)
			{
				Sensor *sNb = *it_nb;
				if (sNb->broadcastCovered == 0 )
				{
					newFrontier->insert(sNb);
					//sNb->broadcastCovered = 1;Never assign this if you want to build a tree.
					touchedSensors.insert(sNb);

					if (!sNb->isChildOnTree)
					{
						sNb->isChildOnTree = true;
						sNb->parentOnTree = s;
						s->childrenOnTree.insert(sNb);
					}
				}
			}
		}

		for (auto itn = newFrontier->begin(); itn != newFrontier->end(); ++ itn)
		{
			Sensor* newS = *itn;
			newS->broadcastCovered = 1;
		}

		delete frontier;
		frontier = newFrontier;
	}

	delete frontier;

	// Reinitialize the flags.
	for (auto it_ts = touchedSensors.begin(); it_ts !=touchedSensors.end();++ it_ts)
	{
		Sensor* s =  *it_ts;

		s->broadcastCovered = 0;
		s->isChildOnTree = false;

	}

}

void Sensor::namingSensors()
{
	// this function is called by a root sensor.
	this->name.clear();
	for (int i = 0; i < field->nSensors; ++i)
	{
		Sensor* s = &(field->sensorPool[i]);
		if (this != s)
		{
			s->name.clear();
			s->name.push_back(s);
			Sensor* sParent = s->parentOnTree;
			while (sParent != this)
			{
				s->name.push_back(sParent);
				sParent = sParent->parentOnTree;
			}
			s->name.push_back(this);
		}
	}
}
int Sensor::distanceOnTree(Sensor* destSensor)
{
	int treeDistance = 0;
	bool jumpOutLoop = false;
	vector<Sensor*> tmpPath;
	if (this == destSensor)
	{
		treeDistance = 0;
	}
	else
	{
		if (this->name.empty())
		{
			treeDistance = destSensor->name.size()-1;
		}
		else if (destSensor->name.empty())
		{
			treeDistance = this->name.size()-1;
		}
		else
		{
			for (int i = 0; i < this->name.size(); ++i)
			{
				for (int j = 0; j < destSensor->name.size();++j)
				{
					if (this->name[i] == destSensor->name[j])
					{
						treeDistance = i + j;
						jumpOutLoop = true;
						break;
					}
				}

				if (jumpOutLoop)
				{
					break;
				}
			}
		}
	}
	assert(treeDistance>=0);
	return treeDistance;
}
Sensor* Sensor::treeNexthop(Sensor* destSensor)
{
	int minDist = 0;
	minDist = distanceOnTree(destSensor);
	Sensor* nextHop = NULL;
	bool hasFind = false;


	for (auto it_nb = neighbors.begin(); it_nb != neighbors.end(); ++ it_nb)
	{
		Sensor* iNb = *it_nb;

		int dist = 0;

		dist = iNb->distanceOnTree(destSensor);

		if (dist < minDist)
		{
			hasFind = true;
			minDist = dist;
			nextHop = iNb;
		}
	}

	if (!hasFind)
	{
		cout<<"There is something wrong with the tree routing..."<<endl;
		cout<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}

	return nextHop;
}
void Sensor::treeRoutingToDest(Sensor* destSensor, vector<Sensor*>& treePath)
{
	Sensor* curSensor = this;
	Sensor* next = NULL;
	while (curSensor != destSensor)
	{
		next = curSensor->treeNexthop(destSensor);
		treePath.push_back(next);
		curSensor = next;
	}
}
void Sensor::sphericalRouting(Sensor* destSensor, vector<Sensor*>& routingPath)
{
	//trivial case
	if (destSensor == this)
		return;
	if (destSensor == NULL)
	{
		cerr<<"Error, the dest sensor should not be null."<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	// Routing on the same type surface.
	int innermostCommmonType = inCommonType(destSensor);
	/*cerr<<"common level "<<innermostCommmonType<<endl;*/
	if (innermostCommmonType != -1)
	{
		sphericalRoutingOnSameType(destSensor,routingPath);
	}
	else
	{
		int srcLowestLevel  = *(this->levelsBelongTo.begin());
		int destLowestLevel = *(destSensor->levelsBelongTo.begin());

// 		cerr<<"Src Lowest level: "<<srcLowestLevel<<" Dest lowest level "<<destLowestLevel<<endl;
// 		cerr<<"Src name size: "<<name.size()<<" dest name size: "<<destSensor->name.size()<<endl;

		if (srcLowestLevel == 0 || destLowestLevel == 0)
		{
			floodingRouteOnFullTopology(destSensor,routingPath);// changed here to avoid tree operation.
		}
		else
		{
			// neither sensors are in the 0 level.
			//first routing on the tree.
			if (srcLowestLevel ==  destLowestLevel)
			{
				cerr<<"The src lowest level should not equal to dest lowest level"<<endl;
				cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
				exit(EXIT_FAILURE);
			}
			
			else
			{
				Sensor* srcLandmark  = field->sensorPool + this->iType2LandmarkIndex[srcLowestLevel];
				Sensor* destLandmark = field->sensorPool + destSensor->iType2LandmarkIndex[destLowestLevel];

				Sensor* newDest = destLandmark;//*(destLandmark->itypeEdgeAttachedLandmarkNb[destLowestLevel].begin());
				//first seg path

				vector<Sensor*>src2Lm;
				this->floodingRouteOnFullTopology(srcLandmark,src2Lm);
				//the second segmentPath;landmark to landmark
				vector<Sensor*>Lm2LmSeg;
				srcLandmark->floodingRouteOnFullTopology(newDest,Lm2LmSeg);

				vector<Sensor*>sameLevelPath;
				newDest->sphericalRoutingOnSameType(destSensor,sameLevelPath);
				/*cout<<"SamelevelPath size: "<<sameLevelPath.size()<<endl;*/
				routingPath.insert(routingPath.end(), src2Lm.begin(),src2Lm.end());
				routingPath.insert(routingPath.end(), Lm2LmSeg.begin(),Lm2LmSeg.end());
				routingPath.insert(routingPath.end(),sameLevelPath.begin(),sameLevelPath.end());


			}
//			if (srcLowestLevel > destLowestLevel)
//			{	
//				Sensor* srcLandmark  = field->sensorPool + iType2LandmarkIndex[srcLowestLevel];
//				Sensor* destLandmark = field->sensorPool + iType2LandmarkIndex[destLowestLevel];
//
//				Sensor* newDest = *(destLandmark->itypeEdgeAttachedLandmarkNb[destLowestLevel].begin());
//				//first seg path
//
//				vector<Sensor*>src2Lm;
//				this->floodingRouteOnFullTopology(srcLandmark,src2Lm);
//				//the second segmentPath;landmark to landmark
//				vector<Sensor*>Lm2LmSeg;
//				srcLandmark->floodingRouteOnFullTopology(newDest,Lm2LmSeg);
//// 
//// 				int  currentLevel = srcLowestLevel;
//// 				
//// 				while(currentLevel > destLowestLevel)
//// 				{
//// 					// 
//// 				}
//// 
//// 
//// 
//// 
//// 				//destLowestLevel  += field->CORE_RADIUS;
//// 				int diffLevel = srcLowestLevel - destLowestLevel;
//// 				vector<Sensor*>treePath;
//// 				// up to the root to route to the destlowestlevel;
//// 				Sensor* sameBranchLevelsensorSrc = name[diffLevel];
//// 				this->treeRoutingToDest(sameBranchLevelsensorSrc,treePath);
//// 				// then routing on the same level
//// 				/*cout<<"treepath size: "<<treePath.size()<<endl;*/
//				vector<Sensor*>sameLevelPath;
//				newDest->sphericalRoutingOnSameType(destSensor,sameLevelPath);
//				/*cout<<"SamelevelPath size: "<<sameLevelPath.size()<<endl;*/
//				routingPath.insert(routingPath.end(), src2Lm.begin(),src2Lm.end());
//				routingPath.insert(routingPath.end(), Lm2LmSeg.begin(),Lm2LmSeg.end());
//				routingPath.insert(routingPath.end(),sameLevelPath.begin(),sameLevelPath.end());
//
//			}
//			else if (srcLowestLevel < destLowestLevel)
//			{
//
//				Sensor* srcLandmark  = field->sensorPool + iType2LandmarkIndex[srcLowestLevel];
//				Sensor* destLandmark = field->sensorPool + iType2LandmarkIndex[destLowestLevel];
//
//				Sensor* newDest = *(destLandmark->itypeEdgeAttachedLandmarkNb[destLowestLevel].begin());
//				//first seg path
//
//				vector<Sensor*>src2Lm;
//				this->floodingRouteOnFullTopology(srcLandmark,src2Lm);
//				//the second segmentPath;landmark to landmark
//				vector<Sensor*>Lm2LmSeg;
//				srcLandmark->floodingRouteOnFullTopology(newDest,Lm2LmSeg);
//				// 
//				// 				int  currentLevel = srcLowestLevel;
//				// 				
//				// 				while(currentLevel > destLowestLevel)
//				// 				{
//				// 					// 
//				// 				}
//				// 
//				// 
//				// 
//				// 
//				// 				//destLowestLevel  += field->CORE_RADIUS;
//				// 				int diffLevel = srcLowestLevel - destLowestLevel;
//				// 				vector<Sensor*>treePath;
//				// 				// up to the root to route to the destlowestlevel;
//				// 				Sensor* sameBranchLevelsensorSrc = name[diffLevel];
//				// 				this->treeRoutingToDest(sameBranchLevelsensorSrc,treePath);
//				// 				// then routing on the same level
//				// 				/*cout<<"treepath size: "<<treePath.size()<<endl;*/
//				vector<Sensor*>sameLevelPath;
//				newDest->sphericalRoutingOnSameType(destSensor,sameLevelPath);
//				/*cout<<"SamelevelPath size: "<<sameLevelPath.size()<<endl;*/
//				routingPath.insert(routingPath.end(), src2Lm.begin(),src2Lm.end());
//				routingPath.insert(routingPath.end(), Lm2LmSeg.begin(),Lm2LmSeg.end());
//				routingPath.insert(routingPath.end(),sameLevelPath.begin(),sameLevelPath.end());
//
//
//				// down to find the child which belong the same level with dest sensor.
//				// down to the child.// flooding to find the dest level sensor.
////				vector<Sensor*> toChildTotalPath;
////				toChildTotalPath.clear();
////				int  currentLevel = srcLowestLevel;
////				Sensor* curS = this;
////				Sensor* next = NULL;
////				while(currentLevel < destLowestLevel)
////				{
////					// down to first child
////					next = curS->toChildPath.back();
////					toChildTotalPath.insert(toChildTotalPath.end(), curS->toChildPath.begin(),curS->toChildPath.end());
////					curS = next;
////					++ currentLevel;
////				}
////				//
////				//cout<<"run wrong here.."<<endl;
////				// cout<< currentLevel<<" "<<destLowestLevel<<endl;
////// 				// srcLowestLevel += field->CORE_RADIUS;
////// 				int diffLevel = destLowestLevel - srcLowestLevel;
////// 				vector<Sensor*>treePath;
////// 				// up to the root to route to the destlowestlevel;
////// 				Sensor* sameBranchLevelsensorDest = destSensor->name[diffLevel];
////				
////				// then routing on the same level
////				vector<Sensor*>sameLevelPath;
////				curS->sphericalRoutingOnSameType(destSensor,sameLevelPath);
////				//sameBranchLevelsensorDest->treeRoutingToDest(destSensor,treePath);
////				routingPath.insert(routingPath.end(),toChildTotalPath.begin(),toChildTotalPath.end());
////				routingPath.insert(routingPath.end(),sameLevelPath.begin(),sameLevelPath.end());
////				
//			}
//
			

		}

	}

}
void Sensor::sphericalRoutingOnSameType(Sensor* destSensor,vector<Sensor*>&sameTypePath)
{
	//trivial case
	if (destSensor == this)
		return;
	if (destSensor == NULL)
	{
		cerr<<"Error, the dest sensor should not be null."<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	// Routing on the same type surface.
	int innermostCommmonType = inCommonType(destSensor);
	if (innermostCommmonType == -1)
	{
		cerr<<"Error, the two sensors are not in the same level."<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}

	if (innermostCommmonType == 0)
	{
		//cout<<"innermost common type..."<<endl;
		floodingRouteOnFullTopology(destSensor,sameTypePath);
	}
	else
	{
		int destLmIndex = destSensor->iType2LandmarkIndex[innermostCommmonType];
		Sensor* destLm = &(field->sensorPool[destLmIndex]);
		int myLmIndex = iType2LandmarkIndex[innermostCommmonType];
		Sensor* myLm = &(field->sensorPool[myLmIndex]);
		if (iType2LandmarkIndex[innermostCommmonType] == destSensor->iType2LandmarkIndex[innermostCommmonType])
		{
			// cout<<"same landmark index..."<<endl;
			floodingRouteOnFullTopology(destSensor,sameTypePath);
		}
		else if (myLm->itypeEdgeAttachedLandmarkNb[innermostCommmonType].find(destLm) != myLm->itypeEdgeAttachedLandmarkNb[innermostCommmonType].end())
		{
			floodingRouteOnFullTopology(destSensor,sameTypePath);
		}
		else
		{

			//Sensor* srcLandmark  = field->sensorPool + this->iType2LandmarkIndex[srcLowestLevel];
			//Sensor* destLandmark = field->sensorPool + destSensor->iType2LandmarkIndex[destLowestLevel];

			//Sensor* newDest = *(destLm->itypeEdgeAttachedLandmarkNb[innermostCommmonType].begin());
			Sensor* newDest = destLm;
			//first seg path

			vector<Sensor*>src2Lm;
			this->floodingRouteOnFullTopology(myLm,src2Lm);
			//the second segmentPath;landmark to landmark
			vector<Sensor*>Lm2LmSeg;
			myLm->floodingRouteOnFullTopology(newDest,Lm2LmSeg);

// 			vector<Sensor*>newDest2DestLm;
// 			newDest->floodingRouteOnFullTopology(destLm,newDest2DestLm);

			vector<Sensor*>lm2Dest;
			destLm->floodingRouteOnFullTopology(destSensor,lm2Dest);
			//vector<Sensor*>sameLevelPath;
			//newDest->sphericalRoutingOnSameType(destSensor,sameLevelPath);
			/*cout<<"SamelevelPath size: "<<sameLevelPath.size()<<endl;*/
// 			routingPath.insert(routingPath.end(), src2Lm.begin(),src2Lm.end());
// 			routingPath.insert(routingPath.end(), Lm2LmSeg.begin(),Lm2LmSeg.end());

			sameTypePath.insert(sameTypePath.end(),src2Lm.begin(),src2Lm.end());

			sameTypePath.insert(sameTypePath.end(),Lm2LmSeg.begin(),Lm2LmSeg.end());
			//sameTypePath.insert(sameTypePath.end(),newDest2DestLm.begin(),newDest2DestLm.end());
			sameTypePath.insert(sameTypePath.end(),lm2Dest.begin(),lm2Dest.end());

			//routingPath.insert(routingPath.end(),sameLevelPath.begin(),sameLevelPath.end());
//			// first route to my landmark.
//// 			int myLmIndex = iType2LandmarkIndex[innermostCommmonType];
//// 			Sensor* myLm = &(field->sensorPool[myLmIndex]);
//
//			//Sensor* ;
//// 			if (myLm == this)
//// 			{
//// 				cout<<"this sensor is a landmark."<<endl;
//// 			}
//			vector<Sensor*>firstSegPath;
//			
//			floodingOnItypeBoundaryRouteTo(innermostCommmonType,myLm,firstSegPath);
//			// cout<<firstSegPath.size()<<" firstseg.size"<<endl;
//			Sensor* lastStep = NULL;
//
//			if (firstSegPath.empty())
//			{
//				//this sensor is a landmark in inner most common type.
//				lastStep = myLm;
//			}
//			else
//			{
//				lastStep = firstSegPath.back();
//			}
//			
//			// then my landmark route to a sensor with the same dest landmark index.
//			// in geographic greedy way.
//			// not proper here for a non-uniform sphere.
//
//			// flooding the the dest landmark on the attached landmark neighbors.
//// 			int destLmIndex = destSensor->iType2LandmarkIndex[innermostCommmonType];
//// 			Sensor* destLm = &(field->sensorPool[destLmIndex]);
//			vector<Sensor*>LMpathLM2LM;
//			lastStep->floodOnItypeEdgeAttachedLandmarkNbsRoutes(innermostCommmonType,destLm,LMpathLM2LM);
//			/*cout<<LMpathLM2LM.size()<<" LMpathLM2LM.size"<<endl;*/
//			vector<Sensor*>actualPathLM2LM;
//			Sensor* formerLM = myLm;
//			for (auto it_iLM = LMpathLM2LM.begin();it_iLM != LMpathLM2LM.end(); ++ it_iLM)
//			{
//				vector<Sensor*>pathNeighboringlandmark;
//				Sensor* latterLM = *it_iLM;
//				formerLM->floodingOnItypeBoundaryRouteTo(innermostCommmonType,latterLM,pathNeighboringlandmark);
//				actualPathLM2LM.insert(actualPathLM2LM.end(),pathNeighboringlandmark.begin(),pathNeighboringlandmark.end());
//				formerLM = latterLM;
//			}
//			/*cout<<actualPathLM2LM.size()<<" actualPathLM2LM.size"<<endl;*/
//			// find the first node that landmark index is equal to the dest landmark index.
//			auto it_firstDestLMIdx = find_if(actualPathLM2LM.begin(),actualPathLM2LM.end(),[&](Sensor* s)->bool
//			{
//				return (s->iType2LandmarkIndex[innermostCommmonType] == destLmIndex);
//			});
//
//			// the last segment.// from the it_first to the destSensor, flood on the itype boundary.
//			Sensor* firstDestLMIdx = *it_firstDestLMIdx;
//			vector<Sensor*> lastSegPath;
//			firstDestLMIdx->floodingOnItypeBoundaryRouteTo(innermostCommmonType,destSensor,lastSegPath);
//			/*cout<<lastSegPath.size()<<" lastSegPath.size"<<endl;*/
//			sameTypePath.insert(sameTypePath.end(),firstSegPath.begin(),firstSegPath.end());
//			sameTypePath.insert(sameTypePath.end(),actualPathLM2LM.begin(),it_firstDestLMIdx + 1);// the last element 
//			sameTypePath.insert(sameTypePath.end(),lastSegPath.begin(),lastSegPath.end());
//

		}

	}
	

}
void Sensor::floodingOnItypeBoundaryRouteTo(const int itype, Sensor* destSensor, vector<Sensor*>& boundPath)
{
	// destSensor should not be null.
	assert(destSensor);
	// This node should be in the same boundary with destSensor.
	assert(isBelongToType(itype) && destSensor->isBelongToType(itype));
	// the boundary path should be empty.
	assert(boundPath.empty());

	int i;
	int hopCount = 0;
	set<int> touchedNodes;

	boundPath.clear();

	//trivial case
	if (destSensor == this)
		return;

	//the path is in the cache
	if (spMapOnItypeBoundary[itype].find(destSensor->index) != spMapOnItypeBoundary[itype].end())
	{
		Sensor* curNode = this;

		while (curNode != destSensor)
		{
			Sensor *curParent = &globalField->sensorPool[curNode->spMapOnItypeBoundary[itype][destSensor->index].second];
			boundPath.push_back(curParent);
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
				for (auto it_sNb = s->neighbors.begin(); it_sNb != s->neighbors.end(); ++ it_sNb)
				{
					Sensor *n = *it_sNb;
					// this condition limited the flooding on the boundary.
					if (n->broadcastCovered == 0 &&  n->isBelongToType(itype))
					{
						n->spMapOnItypeBoundary[itype].insert(DEST2SP::value_type(index, SPLEN_NEXTHOP(hopCount + 1, s->index)));
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
		boundPath.insert(boundPath.begin(), curNode);
		Sensor *curParent = &globalField->sensorPool[curNode->spMapOnItypeBoundary[itype][index].second];
		//also remember the reverse route
		curParent->spMapOnItypeBoundary[itype].insert(DEST2SP::value_type(destSensor->index, SPLEN_NEXTHOP(boundPath.size(), curNode->index)));
		curNode = curParent;
	}
	while (curNode != this);

	for(set<int>::const_iterator iter = touchedNodes.begin(); iter != touchedNodes.end(); ++iter)
	{
		globalField->sensorPool[*iter].broadcastCovered = 0;
		globalField->sensorPool[*iter].parent = NULL;
	}
}
void Sensor::floodingOnBoundaryRouteTo(Sensor* destSensor, vector<Sensor*>& boundPath)
{
	// destSensor should not be null.
	assert(destSensor);
	// This node should be in the same boundary with destSensor.
	assert(boundarySensorType == destSensor->boundarySensorType);
	// the boundary path should be empty.
	assert(boundPath.empty());

	int thisBoundaryType = boundarySensorType;
	int i;
	int hopCount = 0;
	set<int> touchedNodes;

	boundPath.clear();

	//trivial case
	if (destSensor == this)
		return;

	//the path is in the cache
	if (spMapOnBoundary.find(destSensor->index) != spMapOnBoundary.end())
	{
		Sensor* curNode = this;

		while (curNode != destSensor)
		{
			Sensor *curParent = &globalField->sensorPool[curNode->spMapOnBoundary[destSensor->index].second];
			boundPath.push_back(curParent);
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
				for (auto it_sNb = s->neighbors.begin(); it_sNb != s->neighbors.end(); ++ it_sNb)
				{
					Sensor *n = *it_sNb;
					// this condition limited the flooding on the boundary.
					if (n->broadcastCovered == 0 && thisBoundaryType == n->boundarySensorType)
					{
						n->spMapOnBoundary.insert(DEST2SP::value_type(index, SPLEN_NEXTHOP(hopCount + 1, s->index)));
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
		boundPath.insert(boundPath.begin(), curNode);
		Sensor *curParent = &globalField->sensorPool[curNode->spMapOnBoundary[index].second];
		//also remember the reverse route
		curParent->spMapOnBoundary.insert(DEST2SP::value_type(destSensor->index, SPLEN_NEXTHOP(boundPath.size(), curNode->index)));
		curNode = curParent;
	}
	while (curNode != this);

	for(set<int>::const_iterator iter = touchedNodes.begin(); iter != touchedNodes.end(); ++iter)
	{
		globalField->sensorPool[*iter].broadcastCovered = 0;
		globalField->sensorPool[*iter].parent = NULL;
	}
}
void Sensor::floodingOnBoundary(const int itype, const int kHOP,vector<Sensor*>& floodBoundary)
{
	// broadcastCovered are not reinitialized int this function.

	// Flooding should be on the boundary.
	int iTypeBoundary = itype;
	assert(iTypeBoundary > -1);

	// This sensor is a landmark sensor.
	assert(landmark);

	int depth = 0;
	broadcastCovered = 1;

	set<Sensor*> *frontier = new set<Sensor*>;
	frontier->insert(this);

	while (depth != kHOP && ! frontier->empty())
	{
		set<Sensor*> *newFrontier = new set<Sensor*>;
		for (auto it = frontier->begin(); it != frontier->end(); ++ it)
		{
			Sensor* f = *it;
			for (auto fNbit = f->neighbors.begin(); fNbit != f->neighbors.end(); ++ fNbit)
			{
				Sensor* fnb = *fNbit;
				// not broadcast and in the same boundary.
				if (fnb->broadcastCovered == 0 && iTypeBoundary == fnb->boundarySensorType)
				{
					newFrontier->insert(fnb);
				}
			}
		}

		if (!newFrontier->empty())
		{
			++ depth;
		}
		
		for (auto itNF = newFrontier->begin(); itNF != newFrontier->end(); ++ itNF)
		{
			Sensor* sNF = *itNF;
			if (depth != kHOP)
			{	
				sNF->broadcastCovered = 1;
				// sensors in the new frontier are boundary nodes with the same type.
				assert(sNF->boundarySensorType == iTypeBoundary);
				sNF->inactive = true;
				// An active sensor should not be a landmark.
				assert(! sNF->landmark);
			}
			else
			{
				// flooding boundary are not covered, not active.
				floodBoundary.push_back(sNF);
			}
			
		}
		delete frontier;
		frontier = newFrontier;
	}
	delete frontier;

	// Renew the flooding boundary.
	floodBoundary.erase((remove_if(floodBoundary.begin(),floodBoundary.end(),[&](Sensor* lamdaS) -> bool
	{
		return(lamdaS->broadcastCovered == 1);
	})),floodBoundary.end());


}
//for selecting landmarks version2
void Sensor::floodingOnBoundary(const int itype, const int khop)
{
	assert(!iType2Landmark[itype]);
	// this node is set as a landmark.
	iType2Landmark[itype] = true;

	int depth = 0;
	vector<Sensor*>touchedSensors;

	set<Sensor*> *frontier = new set<Sensor*>;
	broadcastCovered = 1;
	touchedSensors.push_back(this);
	
	frontier->insert( this );
	while(depth != khop && !frontier->empty())
	{
		set<Sensor*> *newFrontier = new set<Sensor*>;
		for(set<Sensor*>::const_iterator iFron = frontier->begin(); iFron != frontier->end(); ++iFron) 
		{
			Sensor *f = *iFron;

			for(SCITER NeighborPos = f->neighbors.begin(); NeighborPos != f->neighbors.end(); ++NeighborPos)
			{
				Sensor *s = *NeighborPos;
				if(s->broadcastCovered == 0 && s->isBelongToType(itype))
				{
					newFrontier->insert( s );
				}
			}
		}
		if(!newFrontier->empty())
		{
			++depth;
		}

		//dealing newFrontier
		set<Sensor*>::iterator newFrontierpos;
		for(newFrontierpos = newFrontier->begin(); newFrontierpos != newFrontier->end(); ++newFrontierpos)
		{
			Sensor *s = *newFrontierpos;
			s->broadcastCovered = 1;
			touchedSensors.push_back(s);

			if( s->isBelongToType(itype) && ! s->iType2Landmark[itype])
			{
				s->iType2inactive[itype] = true;
			}

		}
		delete frontier;
		frontier = newFrontier;
	}
	delete frontier;
	
 	int count = touchedSensors.size();
 	if(count < field->avg1HopNeighborhoodSize)
 	{
 		--Number_of_Landmark;
 		this->iType2Landmark[itype] = false;
 		this->iType2inactive[itype] = true;
 	    fprintf_s(stderr,"the touched size is small, cancel this landmark to be an inactive sensor!\n");
 	}

	for(int i = 0; i < touchedSensors.size(); ++i)
	{
		(touchedSensors[i])->broadcastCovered = 0;
	}

}
void Sensor::LocalFlooding(int khop)
{
	int depth = 0;
	set<Sensor*> *frontier = new set<Sensor*>;
	broadcastCovered = 1;
	frontier->insert( this );
	while(depth != khop && !frontier->empty())
	{
		set<Sensor*> *newFrontier = new set<Sensor*>;
		for(set<Sensor*>::const_iterator iFron = frontier->begin(); iFron != frontier->end(); ++iFron) 
		{
			Sensor *f = *iFron;

			for(SCITER NeighborPos = f->neighbors.begin(); NeighborPos != f->neighbors.end(); ++NeighborPos)
			{
				Sensor *s = *NeighborPos;
				if(s->broadcastCovered == 0 && s->onEdge_Degree!=0)
				{
					newFrontier->insert( s );
					//s->broadcastCovered = 1;
				}
			}
		}
		if(!newFrontier->empty())
		{
			++depth;
		}
		
		//dealing newFrontier
		set<Sensor*>::iterator newFrontierpos;
		for(newFrontierpos = newFrontier->begin(); newFrontierpos != newFrontier->end(); ++newFrontierpos)
		{
			Sensor *s = *newFrontierpos;
			s->broadcastCovered = 1;
			if(s->onEdge_Degree != 0 && !s->landmark)
			{
				s->inactive = true;
				if(s->landmark)
				{
					cerr<<"Sensor "<<s->index<<" is a landmark and inactive!!!\n";
					fprintf_s(stderr," Error in %s %d\n",  __FILE__, __LINE__);
					exit(0);
				}
			}
			
		}
		delete frontier;
		frontier = newFrontier;
	}
	delete frontier;
	//Broadcast sensors
	int count = 0;
	for(unsigned int i = 0; i < field->BoundarySensorSet.size(); ++i)
	{
		if(field->BoundarySensorSet[i]->broadcastCovered == 1)
		{
			++count;
		}
	}
	printf("BroadcastCovered sensors %d\n", count);
	if(count < field->avg1HopNeighborhoodSize)
	{
		--Number_of_Landmark;
		this->landmark = false;
		this->inactive = true;
		printf("Cancel this landmark to be an inactive sensor!\n");
	}
	for(int i = 0; i < field->nSensors; ++i)
	{
		field->sensorPool[i].broadcastCovered = 0;
	}
}

void Sensor::LocalFloodingforSkeleton(int khop)
{
	int depth = 0;
	set<Sensor*> *frontier = new set<Sensor*>;
	broadcastCovered = 1;
	frontier->insert( this );
	while(depth != khop && !frontier->empty())
	{
		set<Sensor*> *newFrontier = new set<Sensor*>;
		for(set<Sensor*>::const_iterator iFron = frontier->begin(); iFron != frontier->end(); ++iFron) 
		{
			Sensor *f = *iFron;

			for(SCITER NeighborPos = f->neighbors.begin(); NeighborPos != f->neighbors.end(); ++NeighborPos)
			{
				Sensor *s = *NeighborPos;
				if(s->broadcastCovered == 0/* && s->asSkeleton*/)
				{
					newFrontier->insert( s );
				}
			}
		}
		if(!newFrontier->empty())
		{
			++depth;
		}

		//dealing newFrontier
		set<Sensor*>::iterator newFrontierpos;
		for(newFrontierpos = newFrontier->begin(); newFrontierpos != newFrontier->end(); ++newFrontierpos)
		{
			Sensor *s = *newFrontierpos;
			s->broadcastCovered = 1;
			if(s->onEdge_Degree == 0 && !s->landmark)
			{
				s->inactive = true;
				if(s->landmark)
				{
					printf("Sensor %d is a landmark and inactive!!!\n", s->index);
					exit(0);
				}
			}

		}
		delete frontier;
		frontier = newFrontier;
	}
	delete frontier;
	//Broadcast sensors
	int count = 0;
	for(unsigned int i = 0; i < field->SkeletonNodeSet.size(); ++i)
	{
		if(field->SkeletonNodeSet[i]->broadcastCovered == 1)
		{
			++count;
		}
	}
	printf("BroadcastCovered sensors %d\n", count);
	
	for(int i = 0; i < field->nSensors; ++i)
	{
		field->sensorPool[i].broadcastCovered = 0;
	}
}


int Sensor::LocalFloodingForRefineBoundaryPoints(int khop)
{
	int within_number = 1;
	int depth = 0;
	set<Sensor*> *frontier = new set<Sensor*>;
	broadcastCovered = 1;
	frontier->insert( this );
	while(depth != khop && !frontier->empty())
	{
		set<Sensor*> *newFrontier = new set<Sensor*>;
		for(set<Sensor*>::const_iterator iFron = frontier->begin(); iFron != frontier->end(); ++iFron) 
		{
			Sensor *f = *iFron;

			for(SCITER NeighborPos = f->neighbors.begin(); NeighborPos != f->neighbors.end(); ++NeighborPos)
			{
				Sensor *s = *NeighborPos;
				if(s->broadcastCovered == 0)
				{
					newFrontier->insert( s );
				}
			}
		}
		if(!newFrontier->empty())
		{
			++depth;
		}

		//dealing newFrontier
		set<Sensor*>::iterator newFrontierpos;
		for(newFrontierpos = newFrontier->begin(); newFrontierpos != newFrontier->end(); ++newFrontierpos)
		{
			Sensor *s = *newFrontierpos;

			s->broadcastCovered = 1;
			if(s->onEdge_Degree != 0)
			{
				++within_number;
			}
		}
		delete frontier;
		frontier = newFrontier;
	}
	delete frontier;

	for(int i = 0; i < field->nSensors; ++i)
	{
		field->sensorPool[i].broadcastCovered = 0;
	}
	return within_number;
}

int Sensor::LocalFloodingForRefineConvexPoints(int khop)
{
	int within_number = 1;
	int depth = 0;
	set<Sensor*> *frontier = new set<Sensor*>;
	broadcastCovered = 1;
	frontier->insert( this );
	while(depth != khop && !frontier->empty())
	{
		set<Sensor*> *newFrontier = new set<Sensor*>;
		for(set<Sensor*>::const_iterator iFron = frontier->begin(); iFron != frontier->end(); ++iFron) 
		{
			Sensor *f = *iFron;

			for(SCITER NeighborPos = f->neighbors.begin(); NeighborPos != f->neighbors.end(); ++NeighborPos)
			{
				Sensor *s = *NeighborPos;
				if(s->broadcastCovered == 0)
				{
					newFrontier->insert( s );
				}
			}
		}
		if(!newFrontier->empty())
		{
			++depth;
		}

		//dealing newFrontier
		set<Sensor*>::iterator newFrontierpos;
		for(newFrontierpos = newFrontier->begin(); newFrontierpos != newFrontier->end(); ++newFrontierpos)
		{
			Sensor *s = *newFrontierpos;
			
			s->broadcastCovered = 1;
			if(s->convex)
			{
				++within_number;
			}
		}
		delete frontier;
		frontier = newFrontier;
	}
	delete frontier;
	
	for(int i = 0; i < field->nSensors; ++i)
	{
		field->sensorPool[i].broadcastCovered = 0;
	}
	return within_number;
}

int Sensor::LocalFloodingForRefineConcavePoints(int khop)
{
	int within_number = 1;
	int depth = 0;
	set<Sensor*> *frontier = new set<Sensor*>;
	broadcastCovered = 1;
	frontier->insert( this );
	while(depth != khop && !frontier->empty())
	{
		set<Sensor*> *newFrontier = new set<Sensor*>;
		for(set<Sensor*>::const_iterator iFron = frontier->begin(); iFron != frontier->end(); ++iFron) 
		{
			Sensor *f = *iFron;

			for(SCITER NeighborPos = f->neighbors.begin(); NeighborPos != f->neighbors.end(); ++NeighborPos)
			{
				Sensor *s = *NeighborPos;
				if(s->broadcastCovered == 0)
				{
					newFrontier->insert( s );
				}
			}
		}
		if(!newFrontier->empty())
		{
			++depth;
		}

		//dealing newFrontier
		set<Sensor*>::iterator newFrontierpos;
		for(newFrontierpos = newFrontier->begin(); newFrontierpos != newFrontier->end(); ++newFrontierpos)
		{
			Sensor *s = *newFrontierpos;

			s->broadcastCovered = 1;
			if(s->concave)
			{
				++within_number;
			}
		}
		delete frontier;
		frontier = newFrontier;
	}
	delete frontier;

	for(int i = 0; i < field->nSensors; ++i)
	{
		field->sensorPool[i].broadcastCovered = 0;
	}
	return within_number;
}
void Sensor::constructVoroniTile(const int itype, int offset)
{
	// this node is not a sensor node.
	assert(!this->iType2Landmark[itype]);
	//int thisBoundaryType = itype;
	bool Stop = false;
	bool start_offset = false;
	int OFFSET = 0;
	set<Sensor*>broadcastCoveredSet;
	multimap<double, Sensor*, greater<double> > *frontier = new multimap<double, Sensor*, greater<double> >;

	broadcastCovered = 1;
	broadcastCoveredSet.insert(this);

	frontier->insert( make_pair(this->Rou_p, this) );
	do
	{
		multimap<double, Sensor*, greater<double> > *newFrontier = new multimap<double, Sensor*, greater<double> >;
		for(auto iFron = frontier->begin(); iFron != frontier->end(); ++iFron) 
		{
			Sensor *f = iFron->second;

			for(SCITER NeighborPos = f->neighbors.begin(); NeighborPos != f->neighbors.end(); ++NeighborPos)
			{
				Sensor *s = *NeighborPos;
				if(s->broadcastCovered == 0 && s->isBelongToType(itype))
				{
					newFrontier->insert( make_pair(s->Rou_p, s) );
				}
			}
		}
		if(start_offset)
		{
			++OFFSET;
		}
		//dealing newFrontier
		multimap<double, Sensor*, greater<double> >::iterator newFrontierpos;
		for(newFrontierpos = newFrontier->begin(); newFrontierpos != newFrontier->end(); ++newFrontierpos)
		{
			Sensor *s = newFrontierpos->second;
			if(s->iType2Landmark[itype])
			{	
				this->iType2VoronoiCellID[itype].insert(s->index);
				s->broadcastCovered = 1;
				broadcastCoveredSet.insert(s);
				start_offset = true;
				Stop = true;
			}
			else
			{
				s->broadcastCovered = 1;
				broadcastCoveredSet.insert(s);
			}
		}
		delete frontier;
		frontier = newFrontier;

	}while( (OFFSET != offset || !Stop) && !frontier->empty());
	delete frontier;

	for(auto it = broadcastCoveredSet.begin(); it != broadcastCoveredSet.end(); ++ it)
	{
		Sensor* s = *it;
		s->broadcastCovered = 0;
	}
}
void Sensor::LocalFloodingForEstablishCell(int offset)
{
	bool Stop = false;
	bool start_offset = false;
	int OFFSET = 0;
	multimap<double, Sensor*, greater<double> > *frontier = new multimap<double, Sensor*, greater<double> >;
	broadcastCovered = 1;
	frontier->insert( make_pair(this->onEdge_Degree, this) );
	do
	{
		multimap<double, Sensor*, greater<double> > *newFrontier = new multimap<double, Sensor*, greater<double> >;
		for(multimap<double, Sensor*, greater<double> >::const_iterator iFron = frontier->begin(); iFron != frontier->end(); ++iFron) 
		{
			Sensor *f = iFron->second;

			for(SCITER NeighborPos = f->neighbors.begin(); NeighborPos != f->neighbors.end(); ++NeighborPos)
			{
				Sensor *s = *NeighborPos;
				if(s->broadcastCovered == 0 && s->onEdge_Degree != 0)
				{
					newFrontier->insert( make_pair(s->onEdge_Degree, s) );
				}
			}
		}
		if(start_offset)
		{
			++OFFSET;
		}
		//dealing newFrontier
		multimap<double, Sensor*, greater<double> >::iterator newFrontierpos;
		for(newFrontierpos = newFrontier->begin(); newFrontierpos != newFrontier->end(); ++newFrontierpos)
		{
			Sensor *s = newFrontierpos->second;
			if(s->landmark)
			{
				this->VoronoiCellID.insert(s->index);
				s->broadcastCovered = 1;
				start_offset = true;
				Stop = true;
			}
			else
			{
				s->broadcastCovered = 1;
			}
		}
		delete frontier;
		frontier = newFrontier;

	}while( (OFFSET != offset || !Stop) && !frontier->empty());
	delete frontier;
	for(unsigned int i = 0; i < field->BoundarySensorSet.size(); ++i)
	{
		field->BoundarySensorSet[i]->broadcastCovered = 0;
	}
}

void Sensor::LocalFloodingForEstablishCellforSkeleton(int offset)
{
	bool Stop = false;
	bool start_offset = false;
	int OFFSET = 0;
	multimap<double, Sensor*, greater<double> > *frontier = new multimap<double, Sensor*, greater<double> >;
	broadcastCovered = 1;
	frontier->insert( make_pair(this->onEdge_Degree, this) );
	do
	{
		multimap<double, Sensor*, greater<double> > *newFrontier = new multimap<double, Sensor*, greater<double> >;
		for(multimap<double, Sensor*, greater<double> >::const_iterator iFron = frontier->begin(); iFron != frontier->end(); ++iFron) 
		{
			Sensor *f = iFron->second;

			for(SCITER NeighborPos = f->neighbors.begin(); NeighborPos != f->neighbors.end(); ++NeighborPos)
			{
				Sensor *s = *NeighborPos;
				if(s->broadcastCovered == 0)
				{
					newFrontier->insert( make_pair(s->onEdge_Degree, s) );
				}
			}
		}
		if(start_offset)
		{
			++OFFSET;
		}
		//dealing newFrontier
		multimap<double, Sensor*, greater<double> >::iterator newFrontierpos;
		for(newFrontierpos = newFrontier->begin(); newFrontierpos != newFrontier->end(); ++newFrontierpos)
		{
			Sensor *s = newFrontierpos->second;
			if(s->landmark)
			{
				this->VoronoiCellID.insert(s->index);
				s->broadcastCovered = 1;
				start_offset = true;
				Stop = true;
			}
			else
			{
				s->broadcastCovered = 1;
			}
		}
		delete frontier;
		frontier = newFrontier;
	}while( (OFFSET != offset || !Stop) && !frontier->empty());
	delete frontier;
	for(int i = 0; i < field->nSensors; ++i)
	{
		field->sensorPool[i].broadcastCovered = 0;
	}
}

#if 0
void Sensor::LocalFloodingForEstablishCell(int offset)
{
	bool Stop = false;
	bool start_offset = false;
	int OFFSET = 0;
	set<Sensor*> *frontier = new set<Sensor*>;
	broadcastCovered = 1;
	frontier->insert( this );
	do
	{
		set<Sensor*> *newFrontier = new set<Sensor*>;
		for(set<Sensor*>::const_iterator iFron = frontier->begin(); iFron != frontier->end(); ++iFron) 
		{
			Sensor *f = *iFron;

			for(SCITER NeighborPos = f->neighbors.begin(); NeighborPos != f->neighbors.end(); ++NeighborPos)
			{
				Sensor *s = *NeighborPos;
				if(s->broadcastCovered == 0)
				{
					newFrontier->insert( s );
				}
			}
		}
		if(start_offset)
		{
			++OFFSET;
		}
		//dealing newFrontier
		set<Sensor*>::iterator newFrontierpos;
		for(newFrontierpos = newFrontier->begin(); newFrontierpos != newFrontier->end(); ++newFrontierpos)
		{
			Sensor *s = *newFrontierpos;
			if(s->landmark)
			{
				this->VoronoiCellID.insert(s->index);
				s->broadcastCovered = 1;
				start_offset = true;
				Stop = true;
			}
			else
			{
				s->broadcastCovered = 1;
			}
		}
		delete frontier;
		frontier = newFrontier;
		
	}while( (OFFSET != offset || !Stop) && !frontier->empty());
	delete frontier;
	for(int i = 0; i < field->nSensors; ++i)
	{
		field->sensorPool[i].broadcastCovered = 0;
	}
}
#endif

int Sensor::LocalFloodingforGettingOffsetDistance(Sensor *start, Sensor *end)	//flooding from the new 3rdpoint of a line
{
	int i;
	int offset = 0;
	int hopcountstart = 0;
	int hopcountend = 0;
	bool reach_start = false;
	bool reach_end = false;
	for(i = 0; i < field->nSensors; i++) {
		field->sensorPool[i].broadcastCovered = 0;
		field->sensorPool[i].candidateParents.clear();
		field->sensorPool[i].parent = NULL;
		field->sensorPool[i].children.clear();
	}

	int depth = 0;
	broadcastCovered = 1;
	set<Sensor*> *frontier = new set<Sensor*>;
	frontier->insert(this);

	while (!frontier->empty()) {

		set<Sensor*> *newFrontier = new set<Sensor*>;
		for(set<Sensor*>::const_iterator iFron = frontier->begin(); iFron != frontier->end(); ++iFron) {
			Sensor *s = *iFron;

			for(SCITER neighborpos = s->neighbors.begin(); neighborpos != s->neighbors.end(); ++neighborpos) {		
				Sensor *n = *neighborpos;

				if (n->broadcastCovered == 0) {
					newFrontier->insert(n);
					n->candidateParents.insert(s);
				}
			}
		}

		if (!newFrontier->empty()) {
			depth++;
		}

		//let each node in the new level choose a parent
		for(SCITER iNewFron = newFrontier->begin(); iNewFron != newFrontier->end(); ++iNewFron) 
		{

			Sensor* c = *iNewFron;
			c->broadcastCovered = 1;
			if(c->onEdge_Degree != 0)
			{
				if(c->index == start->index)
				{
					hopcountstart = depth;
					reach_start = true;
				}
				if(c->index == end->index)
				{
					hopcountend = depth;
					reach_end = true;
				}
				if(reach_start && reach_end)
				{
					offset = abs(hopcountstart - hopcountend);
					delete newFrontier;
					delete frontier;
					return offset;
				}
			}
		}
		delete frontier;
		frontier = newFrontier;

	}

	delete frontier;

	for(i = 0; i < field->nSensors; i++) {
		field->sensorPool[i].broadcastCovered = 0;
	}
	return offset;
}

bool Sensor::WithinArea(set<Sensor*> area_sensor_set)
{
	for(SCITER ASSpos = area_sensor_set.begin(); ASSpos != area_sensor_set.end(); ++ASSpos)
	{
		Sensor *a = *ASSpos;
		if(this->index == a->index)
		{
			return true;
		}
	}
	return false;
}

vector<Sensor*> Sensor::PathFlood(Sensor *end)
{	
	int i;
	vector<Sensor*> path;
	path.clear();
	for(i = 0; i < field->nSensors; i++) {
		field->sensorPool[i].broadcastCovered = 0;
		field->sensorPool[i].candidateParents.clear();
		field->sensorPool[i].parent = NULL;
		field->sensorPool[i].children.clear();
	}

	int depth = 0;
	broadcastCovered = 1;
	set<Sensor*> *frontier = new set<Sensor*>;
	frontier->insert(this);

	while (!frontier->empty()) {

		set<Sensor*> *newFrontier = new set<Sensor*>;
		for(set<Sensor*>::const_iterator iFron = frontier->begin(); iFron != frontier->end(); ++iFron) {
			Sensor *s = *iFron;

			for(SCITER neighborpos = s->neighbors.begin(); neighborpos != s->neighbors.end(); ++neighborpos) {		
				Sensor *n = *neighborpos;

				if (n->broadcastCovered == 0 && n->onEdge_Degree != 0) {
					newFrontier->insert(n);
					n->candidateParents.insert(s);
				}
			}
		}

		if (!newFrontier->empty()) {
			depth++;
		}

		//let each node in the new level choose a parent
		for(SCITER iNewFron = newFrontier->begin(); iNewFron != newFrontier->end(); ++iNewFron) {

			Sensor* c = *iNewFron;
			c->broadcastCovered = 1;
			if(c->onEdge_Degree != 0)
			{
				Sensor *parentChosen = NULL;

				//choose a random parent from the candidate set
				int randidx = rand() % c->candidateParents.size();
				int count = 0;

				for(SCITER iCan = c->candidateParents.begin(); iCan != c->candidateParents.end(); ++iCan, ++count) {
					Sensor* can = *iCan;

					if (count == randidx) {
						parentChosen = can;
						break;
					}
				}

				c->parent = parentChosen;
				Sensor* curNode = c;

				if(curNode == end)
				{
					path.push_back(c);

					while (curNode != this)
					{
						Sensor *curParent;

						if (curNode == c) {
							curParent = parentChosen;
						}
						else {
							curParent = curNode->parent;
						}

						path.push_back(curParent);
						curNode = curParent;
					}

					for(i = 0; i < field->nSensors; i++) {
						field->sensorPool[i].broadcastCovered = 0;
					}
					delete newFrontier;
					delete frontier;
					return path;
				}
			}
		}
		delete frontier;
		frontier = newFrontier;
		
	}

	delete frontier;

	for(i = 0; i < field->nSensors; i++) {
		field->sensorPool[i].broadcastCovered = 0;
	}
	return path;
}

vector<Sensor*> Sensor::PathFloodforDirectPath(Sensor *end)
{
	int i;
	vector<Sensor*> path;
	path.clear();
	for(i = 0; i < field->nSensors; i++) {
		field->sensorPool[i].broadcastCovered = 0;
		field->sensorPool[i].candidateParents.clear();
		field->sensorPool[i].parent = NULL;
		field->sensorPool[i].children.clear();
	}

	int depth = 0;
	broadcastCovered = 1;
	set<Sensor*> *frontier = new set<Sensor*>;
	frontier->insert(this);

	while (!frontier->empty()) {

		set<Sensor*> *newFrontier = new set<Sensor*>;
		for(set<Sensor*>::const_iterator iFron = frontier->begin(); iFron != frontier->end(); ++iFron) {
			Sensor *s = *iFron;

			for(SCITER neighborpos = s->neighbors.begin(); neighborpos != s->neighbors.end(); ++neighborpos) {		
				Sensor *n = *neighborpos;

				if (n->broadcastCovered == 0) {
					newFrontier->insert(n);
					n->candidateParents.insert(s);
				}
			}
		}

		if (!newFrontier->empty()) {
			depth++;
		}

		//let each node in the new level choose a parent
		for(SCITER iNewFron = newFrontier->begin(); iNewFron != newFrontier->end(); ++iNewFron) 
		{

			Sensor* c = *iNewFron;
			c->broadcastCovered = 1;
			Sensor *parentChosen = NULL;
			//choose a random parent from the candidate set
			int randidx = rand() % c->candidateParents.size();
			int count = 0;

			for(SCITER iCan = c->candidateParents.begin(); iCan != c->candidateParents.end(); ++iCan, ++count) 
			{
				Sensor* can = *iCan;
				if (count == randidx) {
					parentChosen = can;
					break;
				}
			}

			c->parent = parentChosen;
			Sensor* curNode = c;

			if(curNode == end)
			{
				path.push_back(c);
				while (curNode != this)
				{
					Sensor *curParent;
					if (curNode == c) {
						curParent = parentChosen;
					}
					else {
						curParent = curNode->parent;
					}
					path.push_back(curParent);
					curNode = curParent;
				}
				for(i = 0; i < field->nSensors; i++) {
					field->sensorPool[i].broadcastCovered = 0;
				}
				delete newFrontier;
				delete frontier;
				return path;
			}
		}
		delete frontier;
		frontier = newFrontier;

	}

	delete frontier;

	for(i = 0; i < field->nSensors; i++) {
		field->sensorPool[i].broadcastCovered = 0;
	}
	return path;
}

vector<Sensor*> Sensor::PathCriticalityFlooding(Sensor *end, set<Sensor*> area_sensor_set)
{
	int i;
	bool Incircle = false;
	vector<Sensor*> result;
	vector<VSensor> pathset;
	vector<VSensor> temppathset;
	pathset.clear();
	for(i = 0; i < field->nSensors; i++) {
		field->sensorPool[i].broadcastCovered = 0;
		field->sensorPool[i].candidateParents.clear();
		field->sensorPool[i].parent = NULL;
		field->sensorPool[i].children.clear();
	}

	int depth = 0;
	broadcastCovered = 1;
	set<Sensor*> *frontier = new set<Sensor*>;
	frontier->insert(this);

	while (!frontier->empty()) {

		set<Sensor*> *newFrontier = new set<Sensor*>;
		for(set<Sensor*>::const_iterator iFron = frontier->begin(); iFron != frontier->end(); ++iFron) {
			Sensor *s = *iFron;

			for(SCITER neighborpos = s->neighbors.begin(); neighborpos != s->neighbors.end(); ++neighborpos) {		
				Sensor *n = *neighborpos;
				if (n->broadcastCovered == 0 && n->location.distance(s->location) < Gridinterval+0.5/* && n->onEdge_Degree > 0 && Incircle*/) {
					newFrontier->insert(n);
					n->candidateParents.insert(s);
					s->children.insert(n);
				}
			}
		}

		if (!newFrontier->empty()) {
			depth++;
		}

		//let each node in the new level choose a parent
		for(SCITER iNewFron = newFrontier->begin(); iNewFron != newFrontier->end(); ++iNewFron) {

			Sensor* c = *iNewFron;
			
			c->broadcastCovered = 1;
			Incircle = false;
			for(SCITER area_sensor_set_pos = area_sensor_set.begin(); area_sensor_set_pos != area_sensor_set.end(); ++area_sensor_set_pos)
			{
				Sensor *ass = *area_sensor_set_pos;
				if(ass->index == c->index)
				{
					Incircle = true;
					break;
				}
			}
			if(c->onEdge_Degree != 0 && Incircle)
			{
				Sensor* curNode = c;
				//get all possible path
				if(curNode == end)
				{
					VSensor path;
					path.push_back(curNode);
					pathset.push_back(path);
					temppathset.clear();
					do 
					{
						for(unsigned int j = 0 ; j < pathset.size(); ++j)
						{
							VSensor::iterator pathpos;
							pathpos = pathset[j].end();
							--pathpos;
							curNode = *pathpos;
							for(SCITER canParentspos = curNode->candidateParents.begin(); canParentspos != curNode->candidateParents.end(); ++canParentspos)
							{
								Sensor *canPa = *canParentspos;
								VSensor pt;
								pt.clear();
								for(unsigned int k = 0; k < pathset[j].size(); ++k)
								{
									pt.push_back(pathset[j][k]);
								}
								pt.push_back(canPa);
								temppathset.push_back(pt);
							}
						}
						if(!temppathset.empty())
						{
							pathset = temppathset;
							temppathset.clear();
						}
					} while (curNode != this);
					//find a shortest path
					double sum_criticality = 0xffff;
					for(unsigned int m = 0; m < pathset.size(); ++m)
					{
						double sumofc = 0;
						for(unsigned int n = 0; n < pathset[m].size(); ++n)
						{
							sumofc += pathset[m][n]->Criticality;
						}
						if(sum_criticality > sumofc)
						{
							sum_criticality = sumofc;
							result = pathset[m];
						}
					}

					for(i = 0; i < field->nSensors; i++) {
						field->sensorPool[i].broadcastCovered = 0;
					}
					delete newFrontier;
					delete frontier;
					return result;
				}
			}
		}
		delete frontier;
		frontier = newFrontier;
	}
	delete frontier;

	for(i = 0; i < field->nSensors; i++) {
		field->sensorPool[i].broadcastCovered = 0;
	}
	return result;
}

void Sensor::HoleSensorFlooding(set<Sensor*> hole_sensor_set)
{	
	SCITER hsspos;

	/*
	for(i = 0; i < field->nSensors; i++) {
	field->sensorPool[i].broadcastCovered = 0;
	field->sensorPool[i].candidateParents.clear();
	field->sensorPool[i].parent = NULL;
	field->sensorPool[i].children.clear();
	}*/

	for(hsspos = hole_sensor_set.begin(); hsspos != hole_sensor_set.end(); ++hsspos) 
	{
		Sensor *s = *hsspos;
		s->broadcastCovered = 0;
		s->candidateParents.clear();
		s->parent = NULL;
		s->hole_parent = NULL;
		s->children.clear();
	}

	SCITER HoleSetNeighborpos;
	int depth = 0;
	int CoveredNumber = 0;
	broadcastCovered = 1;
	++CoveredNumber;
	set<Sensor*> *frontier = new set<Sensor*>;
	frontier->insert(this);
	while (CoveredNumber != hole_sensor_set.size()) {

		set<Sensor*> *newFrontier = new set<Sensor*>;
		for(set<Sensor*>::const_iterator iFron = frontier->begin(); iFron != frontier->end(); ++iFron) {
			Sensor *s = *iFron;
			if(s->HoleSetNeighbors.size() > 2)
			{
				for(HoleSetNeighborpos = s->HoleSetNeighbors.begin(); HoleSetNeighborpos != s->HoleSetNeighbors.end(); ++HoleSetNeighborpos) {
					Sensor *n = *HoleSetNeighborpos;
					bool ContainInHole = false;
					for(SCITER hssetpos = hole_sensor_set.begin(); hssetpos != hole_sensor_set.end(); ++hssetpos)
					{
						Sensor *ntemp = *hssetpos;
						if(n->index == ntemp->index)
						{
							ContainInHole = true;
							break;
						}
					}
					if (n->broadcastCovered == 0 && ContainInHole) {
						newFrontier->insert(n);
						n->candidateParents.insert(s);
						//n->broadcastCovered = 1;
					}
				}
				if(newFrontier->size() > 2)
				{
					for(SCITER newpos = newFrontier->begin(); newpos != newFrontier->end(); ++newpos)
					{
						Sensor *cn = *newpos;
						if(cn->BelongInHoleID.size() > 1)
						{
							cn->candidateParents.clear();
							newFrontier->erase(newpos);
							break;
						}
					}
				}

			}
			else
			{
				for(HoleSetNeighborpos = s->HoleSetNeighbors.begin(); HoleSetNeighborpos != s->HoleSetNeighbors.end(); ++HoleSetNeighborpos) {
					Sensor *n = *HoleSetNeighborpos;

					if (n->broadcastCovered == 0) {
						newFrontier->insert(n);
						n->candidateParents.insert(s);
						//n->broadcastCovered = 1;
					}
				}
			}
		}

		if (!newFrontier->empty()) {
			depth++;
		}

		//let each node in the new level choose a parent

		for(SCITER iNewFron = newFrontier->begin(); iNewFron != newFrontier->end(); ++iNewFron) {

			Sensor *c = *iNewFron;
			Sensor *parentChosen = NULL;

			//choose a random parent from the candidate set
			int randidx = rand() % c->candidateParents.size();
			int count = 0;

			for(SCITER iCan = c->candidateParents.begin(); iCan != c->candidateParents.end(); ++iCan, ++count) {
				Sensor* can = *iCan;

				if (count == randidx) {
					parentChosen = can;
					break;
				}
			}

			c->level = depth;
			c->broadcastCovered = 1;
			++CoveredNumber;
			c->hole_parent = parentChosen;
			parentChosen->children.insert(c);
		}
		delete frontier;
		frontier = newFrontier;

	}

	delete frontier;
	for(int i = 0; i < field->nSensors; i++) {
		field->sensorPool[i].broadcastCovered = 0;
	}
}

bool Sensor::SearchInKhopNeighbors(Path *original, Path *pa, int hopnum)
{
	int depth = 0;
	int i = 0;
	bool crossstart = false;
	bool crossend = false;
	set<Sensor*> *frontier = new set<Sensor*>;
	broadcastCovered = 1;
	if(pa != NULL)
	{
		frontier->insert( this );
		while(depth != hopnum && !frontier->empty())
		{
			set<Sensor*> *newFrontier = new set<Sensor*>;
			for(set<Sensor*>::const_iterator iFron = frontier->begin(); iFron != frontier->end(); ++iFron) 
			{
				Sensor *f = *iFron;

				for(SCITER NeighborPos = f->neighbors.begin(); NeighborPos != f->neighbors.end(); ++NeighborPos)
				{
					Sensor *s = *NeighborPos;
					if(s->broadcastCovered == 0 /*&& s->onEdge_Degree != 0*/ /*&& s->VoronoiCellID.size() < 2*/)
					{
						newFrontier->insert( s );
					}
				}
			}
			if(!newFrontier->empty())
			{
				++depth;
			}
			//dealing newFrontier
			set<Sensor*>::iterator newFrontierpos;
			for(newFrontierpos = newFrontier->begin(); newFrontierpos != newFrontier->end(); ++newFrontierpos)
			{
				Sensor *s = *newFrontierpos;
				s->broadcastCovered = 1;
				if(depth == hopnum && s->onEdge_Degree != 0)
				{
					for(set<int>::iterator VIDpos = s->VoronoiCellID.begin(); VIDpos != s->VoronoiCellID.end(); ++VIDpos)
					{
						int VID = *VIDpos;
						if(VID == pa->start->index)
						{
							crossstart = true;
						}
						if(VID == pa->end->index)
						{
							crossend = true;
						}
					}
					if(crossstart && crossend)
					{
						for(i = 0; i < field->nSensors; i++) {
							field->sensorPool[i].broadcastCovered = 0;
						}
						delete newFrontier;
						delete frontier;
						return true;		//neighbors lie in different cells
					}
					else
					{
						for(i = 0; i < field->nSensors; i++) {
							field->sensorPool[i].broadcastCovered = 0;
						}
						delete newFrontier;
						delete frontier;
						return false;
					}
				}
				
			}
			delete frontier;
			frontier = newFrontier;
		}
		delete frontier;
	}


	else
	{
		set<int> unwontedVID;
		unwontedVID.clear();
		frontier->insert( this );
		while(depth != hopnum && !frontier->empty())
		{
			set<Sensor*> *newFrontier = new set<Sensor*>;
			for(set<Sensor*>::const_iterator iFron = frontier->begin(); iFron != frontier->end(); ++iFron) 
			{
				Sensor *f = *iFron;

				for(SCITER NeighborPos = f->neighbors.begin(); NeighborPos != f->neighbors.end(); ++NeighborPos)
				{
					Sensor *s = *NeighborPos;
					if(s->broadcastCovered == 0/* && s->onEdge_Degree != 0*/)
					{
						newFrontier->insert( s );
					}
				}
			}
			if(!newFrontier->empty())
			{
				++depth;
			}
			//dealing newFrontier
			set<Sensor*>::iterator newFrontierpos;
			for(newFrontierpos = newFrontier->begin(); newFrontierpos != newFrontier->end(); ++newFrontierpos)
			{
				Sensor *s = *newFrontierpos;
				s->broadcastCovered = 1;
				if(depth == hopnum/* && s->VoronoiCellID.size() < 4*/ && s->onEdge_Degree != 0)
				{
					for(set<int>::iterator VcellIDpos = s->VoronoiCellID.begin(); VcellIDpos != s->VoronoiCellID.end(); ++VcellIDpos)
					{
						int value = *VcellIDpos;
						if(value != original->start->index && value != original->end->index)
						{
							unwontedVID.insert(value);
						}
						
					}
				}
			}
			multimap<int, int>::iterator connectingpos;
			for(connectingpos = field->ConnectingLandmarkIDPairs.begin(); connectingpos != field->ConnectingLandmarkIDPairs.end(); ++connectingpos)
			{
				int id1 = connectingpos->first;
				int id2 = connectingpos->second;
				bool matchid1 = false;
				bool matchid2 = false;
				for(set<int>::iterator unwantedVIDpos = unwontedVID.begin(); unwantedVIDpos != unwontedVID.end(); ++unwantedVIDpos)
				{
					int va = *unwantedVIDpos;
					if(va == id1)
						matchid1 = true;
					if(va == id2)
						matchid2 = true;
				}
				if(matchid1 && matchid2)
				{
					delete newFrontier;
					delete frontier;
					return true;		//neighbors lie in different cells
				}
			}
			delete frontier;
			frontier = newFrontier;
		}
		delete frontier;
	}
	for(i = 0; i < field->nSensors; i++) {
		field->sensorPool[i].broadcastCovered = 0;
	}
	return false;
}

bool Sensor::NotNearEndsOfPath(Path *ph)
{
	SCITER neighborpos;
	//searching in path's neighbors of start
	for(neighborpos = ph->start->neighbors.begin(); neighborpos != ph->start->neighbors.end(); ++neighborpos)
	{
		Sensor *n = *neighborpos;
		if(n->index == index)
		{
			return false;
		}
	}
	//searching in path's neighbors of end
	for(neighborpos = ph->end->neighbors.begin(); neighborpos != ph->end->neighbors.end(); ++neighborpos)
	{
		Sensor *n = *neighborpos;
		if(n->index == index)
		{
			return false;
		}
	}
	return true;
}

int Sensor::FloodingforGetHopsBetweenStartandEnd(Sensor *end)
{
	int depth = 0;
	set<Sensor*> *frontier = new set<Sensor*>;
	broadcastCovered = 1;
	frontier->insert(this);
	while(!frontier->empty())
	{
		set<Sensor*> *newFrontier = new set<Sensor*>;
		for(set<Sensor*>::const_iterator iFron = frontier->begin(); iFron != frontier->end(); ++iFron) 
		{
			Sensor *f = *iFron;

			for(SCITER NeighborPos = f->neighbors.begin(); NeighborPos != f->neighbors.end(); ++NeighborPos)
			{
				Sensor *s = *NeighborPos;
				if(s->broadcastCovered == 0 && s->onEdge_Degree != 0)
				{
					newFrontier->insert(s);
				}
			}
		}
		if(!newFrontier->empty())
		{
			++depth;
		}
		//dealing newFrontier
		set<Sensor*>::iterator newFrontierpos;
		for(newFrontierpos = newFrontier->begin(); newFrontierpos != newFrontier->end(); ++newFrontierpos)
		{
			Sensor *s = *newFrontierpos;
			s->broadcastCovered = 1;
			if(s == end)
			{
				delete newFrontier;
				delete frontier;
				for(int i = 0; i < field->nSensors; ++i)
				{
					field->sensorPool[i].broadcastCovered = 0;
				}
				return depth;
			}
		}
		delete frontier;
		frontier = newFrontier;
	}
	delete frontier;
	for(int i = 0; i < field->nSensors; ++i)
	{
		field->sensorPool[i].broadcastCovered = 0;
	}
	return depth;
}

int Sensor::FloodingforGetHopsBetweenStartandEndforSkeleton(Sensor *end)
{
	int depth = 0;
	set<Sensor*> *frontier = new set<Sensor*>;
	broadcastCovered = 1;
	frontier->insert(this);
	while(!frontier->empty())
	{
		set<Sensor*> *newFrontier = new set<Sensor*>;
		for(set<Sensor*>::const_iterator iFron = frontier->begin(); iFron != frontier->end(); ++iFron) 
		{
			Sensor *f = *iFron;

			for(SCITER NeighborPos = f->neighbors.begin(); NeighborPos != f->neighbors.end(); ++NeighborPos)
			{
				Sensor *s = *NeighborPos;
				if(s->broadcastCovered == 0 && s->onEdge_Degree == 0)
				{
					newFrontier->insert(s);
				}
			}
		}
		if(!newFrontier->empty())
		{
			++depth;
		}
		//dealing newFrontier
		set<Sensor*>::iterator newFrontierpos;
		for(newFrontierpos = newFrontier->begin(); newFrontierpos != newFrontier->end(); ++newFrontierpos)
		{
			Sensor *s = *newFrontierpos;
			s->broadcastCovered = 1;
			if(s == end)
			{
				delete newFrontier;
				delete frontier;
				for(int i = 0; i < field->nSensors; ++i)
				{
					field->sensorPool[i].broadcastCovered = 0;
				}
				return depth;
			}
		}
		delete frontier;
		frontier = newFrontier;
	}
	delete frontier;
	for(int i = 0; i < field->nSensors; ++i)
	{
		field->sensorPool[i].broadcastCovered = 0;
	}
	return depth;
}
int Sensor::distanceBetweenLandmarks(const int itype,int khop)
{
	int depth = 0;
	set<Sensor*> *frontier = new set<Sensor*>;
	vector<Sensor*> touchedSensors;

	broadcastCovered = 1;
	touchedSensors.push_back(this);
	
	frontier->insert( this );

	while(depth != khop && !frontier->empty())
	{
		set<Sensor*> *newFrontier = new set<Sensor*>;
		for(set<Sensor*>::const_iterator iFron = frontier->begin(); iFron != frontier->end(); ++iFron) 
		{
			Sensor *f = *iFron;

			for(SCITER NeighborPos = f->neighbors.begin(); NeighborPos != f->neighbors.end(); ++NeighborPos)
			{
				Sensor *s = *NeighborPos;
				if(s->broadcastCovered == 0 && s->isBelongToType(itype))
				{
					newFrontier->insert( s );
				}
			}
		}
		if(!newFrontier->empty())
		{
			++depth;
		}
		//dealing newFrontier
		set<Sensor*>::iterator newFrontierpos;
		for(newFrontierpos = newFrontier->begin(); newFrontierpos != newFrontier->end(); ++newFrontierpos)
		{
			Sensor *s = *newFrontierpos;

			s->broadcastCovered = 1;
			touchedSensors.push_back(s);

			if(s->iType2Landmark[itype])
			{
				delete newFrontier;
				delete frontier;

				// This step is very important before return.
				for(int i = 0; i < touchedSensors.size(); ++i)
				{
					touchedSensors[i]->broadcastCovered = 0;
				}

				return depth;
			}
		}
		delete frontier;
		frontier = newFrontier;
	}
	delete frontier;

	
	for(int i = 0; i < touchedSensors.size(); ++i)
	{
		touchedSensors[i]->broadcastCovered = 0;
	}

	return depth;

}
int Sensor::LandmarkFlooding(int khop)
{
	int depth = 0;
	set<Sensor*> *frontier = new set<Sensor*>;
	broadcastCovered = 1;
	frontier->insert( this );
	while(depth != khop && !frontier->empty())
	{
		set<Sensor*> *newFrontier = new set<Sensor*>;
		for(set<Sensor*>::const_iterator iFron = frontier->begin(); iFron != frontier->end(); ++iFron) 
		{
			Sensor *f = *iFron;

			for(SCITER NeighborPos = f->neighbors.begin(); NeighborPos != f->neighbors.end(); ++NeighborPos)
			{
				Sensor *s = *NeighborPos;
				if(s->broadcastCovered == 0 && s->onEdge_Degree != 0)
				{
					newFrontier->insert( s );
				}
			}
		}
		if(!newFrontier->empty())
		{
			++depth;
		}
		//dealing newFrontier
		set<Sensor*>::iterator newFrontierpos;
		for(newFrontierpos = newFrontier->begin(); newFrontierpos != newFrontier->end(); ++newFrontierpos)
		{
			Sensor *s = *newFrontierpos;
			s->broadcastCovered = 1;
			if(s->landmark)
			{
				delete newFrontier;
				delete frontier;
				return depth;
			}
		}
		delete frontier;
		frontier = newFrontier;
	}
	delete frontier;
	for(int i = 0; i < field->nSensors; ++i)
	{
		field->sensorPool[i].broadcastCovered = 0;
	}
	return depth;
}

bool Sensor::PathByBoundarySensors(Sensor *end)
{
	int i;
	vector<Sensor*> path;
	path.clear();
	for(i = 0; i < field->nSensors; i++) {
		field->sensorPool[i].broadcastCovered = 0;
		field->sensorPool[i].candidateParents.clear();
		field->sensorPool[i].parent = NULL;
		field->sensorPool[i].children.clear();
	}
	int depth = 0;
	
	set<Sensor*> *frontier = new set<Sensor*>;
	broadcastCovered = 1;
	frontier->insert( this );
	while (!frontier->empty()) {

		set<Sensor*> *newFrontier = new set<Sensor*>;
		for(set<Sensor*>::const_iterator iFron = frontier->begin(); iFron != frontier->end(); ++iFron) {
			Sensor *s = *iFron;

			for(SCITER neighborpos = s->neighbors.begin(); neighborpos != s->neighbors.end(); ++neighborpos) {		
				Sensor *n = *neighborpos;

				if (n->broadcastCovered == 0) {
					newFrontier->insert(n);
					n->candidateParents.insert(s);
				}
			}
		}

		if (!newFrontier->empty()) {
			depth++;
		}

		//let each node in the new level choose a parent
		for(SCITER iNewFron = newFrontier->begin(); iNewFron != newFrontier->end(); ++iNewFron) {

			Sensor* c = *iNewFron;
			c->broadcastCovered = 1;
			Sensor *parentChosen = NULL;

			//choose a random parent from the candidate set
			int randidx = rand() % c->candidateParents.size();
			int count = 0;

			for(SCITER iCan = c->candidateParents.begin(); iCan != c->candidateParents.end(); ++iCan, ++count) {
				Sensor* can = *iCan;
				if (count == randidx) {
					parentChosen = can;
					break;
				}
			}

			c->parent = parentChosen;
			Sensor* curNode = c;

			if(curNode == end)
			{
				path.push_back(c);
				
				while (curNode != this)
				{
					Sensor *curParent;
					if (curNode == c) {
						curParent = parentChosen;
					}
					else {
						curParent = curNode->parent;
					}
					path.push_back(curParent);
					if(curParent->onEdge_Degree == 0)
					{
						return false;
					}
					curNode = curParent;
				}
			}
		}
		delete frontier;
		frontier = newFrontier;
	}

	delete frontier;

	for(i = 0; i < field->nSensors; i++) {
		field->sensorPool[i].broadcastCovered = 0;
	}
	return true;
}

bool Sensor::UnknowDistTo(Sensor *end)
{
	if(this->KnownHopDistEndID.empty())
	{
		printf("KnownHopDistEndID is empty\n");
		return true;
		//exit(0);
	}
	else
	{
		for(unsigned int i = 0; i < this->KnownHopDistEndID.size(); ++i)
		{
			if(end->index == this->KnownHopDistEndID[i])
			{
				return false;
			}
		}
		return true;
	}
}

bool Sensor::FloodingtoFindConcavepointsWithin(int dist, set<Sensor*> &concaveset)
{
	int depth = 0;
	bool result = false;
	set<Sensor*> *frontier = new set<Sensor*>;
	broadcastCovered = 1;
	frontier->insert( this );
	while(depth != dist && !frontier->empty())
	{
		set<Sensor*> *newFrontier = new set<Sensor*>;
		for(set<Sensor*>::const_iterator iFron = frontier->begin(); iFron != frontier->end(); ++iFron) 
		{
			Sensor *f = *iFron;

			for(SCITER NeighborPos = f->neighbors.begin(); NeighborPos != f->neighbors.end(); ++NeighborPos)
			{
				Sensor *s = *NeighborPos;
				if(s->broadcastCovered == 0)
				{
					newFrontier->insert( s );
				}
			}
		}
		if(!newFrontier->empty())
		{
			++depth;
		}

		//dealing newFrontier
		set<Sensor*>::iterator newFrontierpos;
		for(newFrontierpos = newFrontier->begin(); newFrontierpos != newFrontier->end(); ++newFrontierpos)
		{
			Sensor *s = *newFrontierpos;

			s->broadcastCovered = 1;
			if(s->concave)
			{
				result = true;
				concaveset.insert(s);
			}
		}
		delete frontier;
		frontier = newFrontier;
	}
	delete frontier;

	for(int i = 0; i < field->nSensors; ++i)
	{
		field->sensorPool[i].broadcastCovered = 0;
	}
	return result;
}

bool Sensor::BelongOtherIntersectLine(int ID1, int ID2)
{
	multimap<int, int>::iterator endsIDpos;
	if(this->line_ends_ID.empty())
	{
		return false;
	}
	for (endsIDpos = this->line_ends_ID.begin(); endsIDpos != this->line_ends_ID.end(); ++endsIDpos)
	{
		int value1 = endsIDpos->first;
		int value2 = endsIDpos->second;
		if(value1 != ID1 && value1 != ID2 && value2 != ID1 && value2 != ID2)
		{
			return true;
		}
	}
	return false;
}

bool Sensor::NotConnectWith(Sensor *end)
{
	for(set<int>::iterator CIDpos = this->ConnectedID.begin(); CIDpos != this->ConnectedID.end(); ++CIDpos)
	{
		int CID = *CIDpos;
		if(CID == end->index)
		{
			return false;
		}
	}
	return true;
}

void Sensor::LocalFloodingForMeasureTheHopsToSensors(int khop)
{
	int depth = 0;
	set<Sensor*> *frontier = new set<Sensor*>;
	broadcastCovered = 1;
	frontier->insert( this );
	while(depth != khop && !frontier->empty())
	{
		set<Sensor*> *newFrontier = new set<Sensor*>;
		for(set<Sensor*>::const_iterator iFron = frontier->begin(); iFron != frontier->end(); ++iFron) 
		{
			Sensor *f = *iFron;

			for(SCITER NeighborPos = f->neighbors.begin(); NeighborPos != f->neighbors.end(); ++NeighborPos)
			{
				Sensor *s = *NeighborPos;
				if(s->broadcastCovered == 0)
				{
					newFrontier->insert( s );
				}
			}
		}
		if(!newFrontier->empty())
		{
			++depth;
		}

		//dealing newFrontier
		set<Sensor*>::iterator newFrontierpos;
		for(newFrontierpos = newFrontier->begin(); newFrontierpos != newFrontier->end(); ++newFrontierpos)
		{
			Sensor *s = *newFrontierpos;
			s->broadcastCovered = 1;
			if(s->onEdge_Degree == 0)
			{
				s->DistanceToLandmarks.insert(make_pair(depth, this->index));
			}
		}
		delete frontier;
		frontier = newFrontier;
	}
	delete frontier;
	for(int i = 0; i < field->nSensors; ++i)
	{
		field->sensorPool[i].broadcastCovered = 0;
	}
}

vector<Sensor*> Sensor::SpreadPartID(int part_id)
{
	vector<Sensor*> result;
	int count = 0;
	int depth = 0;
	set<Sensor*> *frontier = new set<Sensor*>;
	broadcastCovered = 1;
	this->boundary_partID = part_id;
	frontier->insert( this );
	while(!frontier->empty())
	{
		set<Sensor*> *newFrontier = new set<Sensor*>;
		for(set<Sensor*>::const_iterator iFron = frontier->begin(); iFron != frontier->end(); ++iFron) 
		{
			Sensor *f = *iFron;
			for(SCITER faSPos = f->LandmarkNeighbors.begin(); faSPos != f->LandmarkNeighbors.end(); ++faSPos)
			{
				Sensor *s = *faSPos;
				if(s->broadcastCovered == 0)
				{
					newFrontier->insert( s );
				}
			}
		}
		if(!newFrontier->empty())
		{
			++depth;
		}
		//dealing newFrontier
		set<Sensor*>::iterator newFrontierpos;
		for(newFrontierpos = newFrontier->begin(); newFrontierpos != newFrontier->end(); ++newFrontierpos)
		{
			Sensor *s = *newFrontierpos;
			s->broadcastCovered = 1;
			s->boundary_partID = part_id;
			result.push_back(s);
			++count;
		}
		delete frontier;
		frontier = newFrontier;
	}
	delete frontier;
	cout << "Number of nodes with part ID equal " << part_id << " is: " << count << endl;
	return result;
}