#include "geometry3D.h"
#include "sensor.h"
#include "common.h"
#include <algorithm>
#include "assert.h"
#include "random.h"



/////////////////////////////// Line segments /////////////////////////
LineSegment::LineSegment(int LineID_, Sensor *start_, Sensor *end_)
{
	usedtimes = 0;
	Processingtimes = 0;
	ProcessingLater = false;
	RefinePhaseCanceled = false;
	The3rdPoint.clear();
	LineID = LineID_;
	start = start_;
	end = end_;

	boundaryPath.clear();
	line_Field = NULL;
	valid = false;
	attachedTriangles.clear();
	broadcastedCovered = 0;
	if (start_!= NULL || end_ !=NULL)
	{
		lineImportance = start_->Rou_p + end_->Rou_p;
	}
	
	flowLength = 0.0;
	atBoundary = false;
	CClockEnd = NULL;
	CClockStart = NULL;

	itypeBoundaryPath.clear();
	itypeEdgeAttachedTriangles.clear();
	iTypeFlowLength.clear();
	
}

LineSegment::LineSegment()
{
	usedtimes = 0;
	Processingtimes = 0;
	ProcessingLater = false;
	RefinePhaseCanceled = false;
	The3rdPoint.clear();
	LineID = 0;
	start = NULL;
	end = NULL;

	boundaryPath.clear();
	line_Field = NULL;
	valid = false;
	attachedTriangles.clear();
	broadcastedCovered = 0;
	lineImportance = 0;
	flowLength = 0.0;
	atBoundary = false;

	CClockEnd = NULL;
	CClockStart = NULL;

	itypeBoundaryPath.clear();
	itypeEdgeAttachedTriangles.clear();

	iTypeFlowLength.clear();
}
LineSegment::LineSegment(const LineSegment& otherLCopy)
{
	usedtimes = otherLCopy.usedtimes;
	Processingtimes = otherLCopy.Processingtimes;
	ProcessingLater = otherLCopy.ProcessingLater;
	RefinePhaseCanceled = otherLCopy.RefinePhaseCanceled;
	The3rdPoint = otherLCopy.The3rdPoint;
	LineID = otherLCopy.LineID;
	start = otherLCopy.start;
	end = otherLCopy.end;

	// the push or insert may use the copy function.
	boundaryPath = otherLCopy.boundaryPath;
	itypeBoundaryPath = (otherLCopy.itypeBoundaryPath) ;

	line_Field = otherLCopy.line_Field;
	valid = otherLCopy.valid;
	broadcastedCovered = otherLCopy.broadcastedCovered;
	lineImportance = otherLCopy.lineImportance;
	flowLength = otherLCopy.flowLength;
	atBoundary = otherLCopy.atBoundary;

	CClockEnd = otherLCopy.CClockEnd;
	CClockStart = otherLCopy.CClockStart;

	iTypeFlowLength = otherLCopy.iTypeFlowLength;

}

LineSegment::~LineSegment()
{
}

double LineSegment::distanceToPoint(Point pt)
{
	double x1 = start->location.x;
	double y1 = start->location.y;
	double x2 = end->location.x;
	double y2 = end->location.y;
	double x0 = pt.x;
	double y0 = pt.y;

	double res = (x2-x1)*(y1-y0) - (x1-x0)*(y2-y1);
	res = myabs(res);

	res /= sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));

	return res;
}

bool LineSegment::containEnd(Point pt)
{
	return (start->location.approxEqual(pt) || end->location.approxEqual(pt));
}

bool LineSegment::TheSameLine(LineSegment *l2)
{
	if((start == l2->start && end == l2->end) || (start == l2->end && end == l2->start))
		return true;
	else
		return false;
}

bool LineSegment::CanFindASensorInThe3rdPoint(Sensor *s)
{
	for(unsigned int i = 0; i < The3rdPoint.size(); ++i)
	{
		if(The3rdPoint[i]->index == s->index)
		{
			printf("It's in the 3rd point\n");
			return true;
		}
	}
	return false;
}
//Never connect the corresponding 3rd points which is that both of them belong to the same line's 3rd point
bool LineSegment::CorrespThe3rdPointOfLineVertex(Sensor *s)
{
	unsigned int i;
	for(i = 0; i < start->Corresp3rdPoint.size(); ++i)
	{
		if(start->Corresp3rdPoint[i]->index == s->index)
		{
			printf("Find coorespond point\n");
			return true;
		}
	}
	for(i = 0; i < end->Corresp3rdPoint.size(); ++i)
	{
		if(end->Corresp3rdPoint[i]->index == s->index)
		{
			printf("Find coorespond point\n");
			return true;
		}
	}
	return false;
}

void LineSegment::init(int LineID_, Sensor *start_, Sensor *end_)
{
	usedtimes = 0;
	Processingtimes = 0;
	ProcessingLater = false;
	RefinePhaseCanceled = false;
	The3rdPoint.clear();
	LineID = LineID_;
	start = start_;
	end = end_;
}
Sensor* LineSegment::sharedVertice(LineSegment& otherL)
{
	// this line is not equal to the otherL.
	assert(! TheSameLine(&otherL));
	Sensor* sharedV = NULL;

	if (start == otherL.start || start == otherL.end )
	{
		sharedV = start;
	}
	else if ( end == otherL.end || end == otherL.start)
	{
		sharedV = end;
	}
	
	return sharedV;
	
}
bool LineSegment::isSharingAVertexWith(LineSegment& otherL)
{
	// this line is not equal to the otherL.
	assert(! TheSameLine(&otherL));
	if (LineID == otherL.LineID)
	{
		cerr<<"Error while calling isSharingAVertexWith."<<endl;
		cerr<<"The two lines should not equal."<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<endl;
		exit(EXIT_FAILURE);
	}
	else
	{
		if (start == otherL.start || start == otherL.end || end == otherL.end || end == otherL.start)
		{
			return true;
		}
		else
			return false;
	}
	

}
void LineSegment::floodStart2EndPathOnItypeBoundary(const int itype)
{
	// not a null path
	assert(start != NULL);
	assert(end != NULL);
	// should be on the same boundary
	assert(start->isBelongToType(itype) && end->isBelongToType(itype));

	start->floodingOnItypeBoundaryRouteTo(itype,end,itypeBoundaryPath[itype]);
}
void LineSegment::floodStart2EndPathOnBoundary()
{
	// not a null path
	assert(start != NULL);
	assert(end != NULL);
	// should be on the same boundary
	assert(start->boundarySensorType == end->boundarySensorType);

	start->floodingOnBoundaryRouteTo(end,boundaryPath);

}
void LineSegment::ProcessingThe3rdPoint()
{
	if(The3rdPoint.size() == 2)
	{
		The3rdPoint[0]->Corresp3rdPoint.push_back(The3rdPoint[1]);
		The3rdPoint[1]->Corresp3rdPoint.push_back(The3rdPoint[0]);
	}
}
//
Sensor* LineSegment::GetTheGreatest3rdPoint(multimap< double, Sensor*, greater<double> > sharepoints)
{
	SCITER The3rdpointNeighborpos;
	Sensor *Greatest = NULL;
	bool ContainInLandmarkNeighborhood = false;
	multimap< double, Sensor*, greater<double> > Greatestsamples;
	multimap< double, Sensor*, greater<double> > Landmarksamples;
	multimap< double, Sensor*, greater<double> > Concavesamples;
	//multimap< int, Sensor*, less<double> > CandidateBestResults;
	//multimap< int, Sensor*, less<double> > CandidateResults;
	multimap< double, Sensor*, greater<double> > CandidateBestResults;
	multimap< double, Sensor*, greater<double> > CandidateResults;

	multimap< double, Sensor*, greater<double> >::iterator POS;

	for(POS = sharepoints.begin(); POS != sharepoints.end(); ++POS)
	{
		Sensor *p = POS->second;
		if(start->PathByBoundarySensors(p) && end->PathByBoundarySensors(p))
		{
			CandidateBestResults.insert(make_pair(p->Rou_p, p));
		}
		else
		{
			CandidateResults.insert(make_pair(p->Rou_p, p));
			printf("Detect path through inside\n");
		}
	}
	if(!CandidateBestResults.empty())
	{
		POS = CandidateBestResults.begin();
		while(POS != CandidateBestResults.end())
		{
			Sensor *s = POS->second;
			if(!CanFindASensorInThe3rdPoint(s))
			{
				printf("Get a Greatest result!\n");
				return s;
			}
			++POS;
		}
		
		return POS->second;
	}
	else if (!CandidateResults.empty())
	{
		POS = CandidateResults.begin();
		while(POS != CandidateResults.end())
		{
			Sensor *s = POS->second;
			if(!CanFindASensorInThe3rdPoint(s))
			{
				printf("Get a common result!\n");
				return s;
			}
			++POS;
		}
		return POS->second;
	}
	else
	{
		printf("Find a greatest point part wrong!\n");
		exit(0);
	}
}

Sensor* LineSegment::GetTheGreatest3rdPointforSkeleton(multimap< double, Sensor*, greater<double> > sharepoints)
{
	SCITER The3rdpointNeighborpos;
	Sensor *Greatest = NULL;
	bool ContainInLandmarkNeighborhood = false;
	multimap< double, Sensor*, greater<double> > Greatestsamples;
	multimap< double, Sensor*, greater<double> > Landmarksamples;
	multimap< double, Sensor*, greater<double> > Concavesamples;
	//multimap< int, Sensor*, less<double> > CandidateBestResults;
	//multimap< int, Sensor*, less<double> > CandidateResults;
	multimap< double, Sensor*, greater<double> > CandidateBestResults;
	multimap< double, Sensor*, greater<double> > CandidateResults;

	multimap< double, Sensor*, greater<double> >::iterator POS;

	for(POS = sharepoints.begin(); POS != sharepoints.end(); ++POS)
	{
		Sensor *p = POS->second;
		if(start->PathByBoundarySensors(p) && end->PathByBoundarySensors(p))
		{
			CandidateBestResults.insert(make_pair(p->Rou_p, p));
		}
		else
		{
			CandidateResults.insert(make_pair(p->Rou_p, p));
			printf("Detect path through inside\n");
		}
	}
	if(!CandidateBestResults.empty())
	{
		POS = CandidateBestResults.begin();
		while(POS != CandidateBestResults.end())
		{
			Sensor *s = POS->second;
			if(!CanFindASensorInThe3rdPoint(s))
			{
				printf("Get a Greatest result!\n");
				return s;
			}
			++POS;
		}

		return POS->second;
	}
	else if (!CandidateResults.empty())
	{
		POS = CandidateResults.begin();
		while(POS != CandidateResults.end())
		{
			Sensor *s = POS->second;
			if(!CanFindASensorInThe3rdPoint(s))
			{
				printf("Get a common result!\n");
				return s;
			}
			++POS;
		}
		return POS->second;
	}
	else
	{
		printf("Find a greatest point part wrong!\n");
		exit(0);
	}
}

//if The3rdPoint set is empty, return the first element, otherwise return the point which is not belong the 3rd point's landmark neighborhood
#if 0
Sensor* LineSegment::GetTheGreatest3rdPoint(multimap< double, Sensor*, greater<double> > sharepoints)
{
	SCITER The3rdpointNeighborpos;
	int offset = 100;
	Sensor *Greatest = NULL;
	bool ContainInLandmarkNeighborhood = false;
	multimap< double, Sensor*, greater<double> > Greatestsamples;
	multimap< double, Sensor*, greater<double> > Landmarksamples;

	multimap< double, Sensor*, greater<double> >::iterator POS;

	POS = sharepoints.begin();
	if(The3rdPoint.empty())
	{
		
		for(POS = sharepoints.begin(); POS != sharepoints.end(); ++POS)
		{
			Sensor *s = POS->second;
			int temp_offset = s->LocalFloodingforGettingOffsetDistance(start, end);
			if(offset > temp_offset)
			{
				offset = temp_offset;
				Greatest = s;
			}
		}
		return Greatest;
	}
	else
	{
		for(unsigned int i = 0; i < The3rdPoint.size(); ++i)
		{
			
			for(POS = sharepoints.begin(); POS != sharepoints.end(); ++POS)
			{
				Sensor *g = POS->second;
				for(The3rdpointNeighborpos = The3rdPoint[i]->LandmarkNeighbors.begin(); The3rdpointNeighborpos != The3rdPoint[i]->LandmarkNeighbors.end();
					++The3rdpointNeighborpos)
				{
					Sensor *sen = *The3rdpointNeighborpos;
					
					if(POS->second->index != sen->index)
						continue;
					else
					{
						ContainInLandmarkNeighborhood = true;
						if(!CanFindASensorInThe3rdPoint(g))
						{
							Landmarksamples.insert(make_pair(POS->second->onEdge_Degree, POS->second));
						}
					}
				}
				if(!ContainInLandmarkNeighborhood && !CanFindASensorInThe3rdPoint(g))
					Greatestsamples.insert(make_pair(POS->second->onEdge_Degree, POS->second));
			}
		}
		if(!Greatestsamples.empty())
		{
			Sensor *gr = NULL;
			if(Greatestsamples.size() == 1)
			{
				POS = Greatestsamples.begin();
				gr = POS->second;
				return gr;
			}
			else
			{
				for(POS = Greatestsamples.begin(); POS != Greatestsamples.end(); ++POS)
				{
					gr = POS->second;
					int temp_offset = gr->LocalFloodingforGettingOffsetDistance(start, end);
					if(offset > temp_offset)
					{
						offset = temp_offset;
						Greatest = gr;
					}
				}
			}
			printf("Get a Greatestsamples\n");
			return Greatest;
		}
		else if(!Landmarksamples.empty())
		{
			Sensor *gr = NULL;
			if(Landmarksamples.size() == 1)
			{
				POS = Landmarksamples.begin();
				gr = POS->second;
				return gr;
			}
			else
			{
				for(POS = Landmarksamples.begin(); POS != Landmarksamples.end(); ++POS)
				{
					gr = POS->second;
					int temp_offset = gr->LocalFloodingforGettingOffsetDistance(start, end);
					if(offset > temp_offset)
					{
						offset = temp_offset;
						Greatest = gr;
					}
				}
			}
			printf("Get a Landmarkneighborsamples\n");
			return Greatest;
		}
		else
			return NULL;
		printf("Wrong!\n");
		exit(0);
	}
}
#endif


//////////////////////////////////////////////////////////////////////////Triangle//////////////////////////////////////////////////////////////////////////
Triangle::Triangle(int _Tri_ID, Sensor *t1_, Sensor *t2_, Sensor *t3_)
{
	Tri_ID = _Tri_ID;
	t1 = t1_;
	t2 = t2_;
	t3 = t3_;
	Tri_field = NULL;
	Normalx = 0;
	Normaly = 0;
	Normalz = 0;

	triImportance = t1_->Rou_p + t2_->Rou_p + t3_->Rou_p;
	broadcastCovered = 0;
	
	v2vMap.clear();
	v2theta.clear();
}
Sensor* Triangle::get3rdSensorByEdge(LineSegment* oppEdge)
{
	// the Edge should be a part of the triangle.
	set<Sensor*>vertex;
	Sensor* oppEdgeS = oppEdge->start;
	Sensor* oppEdgeE = oppEdge->end;

	vertex.insert(t1);
	vertex.insert(t2);
	vertex.insert(t3);

	vertex.erase(oppEdgeS);
	vertex.erase(oppEdgeE);

	if (vertex.size() == 1)
	{
		return (*(vertex).begin());
	}
	else
	{
		cerr<<"The edge is not a part of the triangle!"<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
}
LineSegment* Triangle::findASharingItypeEdgeWith(const int itype, Triangle* ls)
{
	LineSegment* sharingEdge = NULL;
	if (this == ls)
	{
		cerr<<"The triangle is  equal to me...!"<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	Field* f = t1->field;
	set<LineSegment*>bothEdges;
	vector<LineSegment*>myEdges;
	vector<LineSegment*>yourEdges;

	auto l1 = f->findItypeEdgePointer(itype,t1,t2);
	auto l2 = f->findItypeEdgePointer(itype,t2,t3);
	auto l3 = f->findItypeEdgePointer(itype,t1,t3);

	myEdges.push_back(l1);
	myEdges.push_back(l2);
	myEdges.push_back(l3);

	auto s1 = (f->findItypeEdgePointer(itype,ls->t1,ls->t2));
	auto s2 = (f->findItypeEdgePointer(itype,ls->t2,ls->t3));
	auto s3 = (f->findItypeEdgePointer(itype,ls->t1,ls->t3));

	yourEdges.push_back(s1);
	yourEdges.push_back(s2);
	yourEdges.push_back(s3);

	bothEdges.insert(l1);
	bothEdges.insert(l2);
	bothEdges.insert(l3);
	bothEdges.insert(s1);
	bothEdges.insert(s2);
	bothEdges.insert(s3);

	if (bothEdges.size() != 5)
	{
		cerr<<"The triangle T and I do not share an edge.!"<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	else
	{
		for (auto itMy = myEdges.begin(); itMy != myEdges.end(); ++ itMy)
		{
			LineSegment* meE = *itMy;
			auto itf = find(yourEdges.begin(),yourEdges.end(),meE);
			if (itf != yourEdges.end())
			{
				sharingEdge = meE;
				break;
			}
		}
	}

	return sharingEdge;


}
LineSegment* Triangle::findASharingEdgeWith(Triangle* ls)
{
	LineSegment* sharingEdge = NULL;
	if (this == ls)
	{
		cerr<<"The triangle is  equal to me...!"<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	Field* f = t1->field;
	set<LineSegment*>bothEdges;
	vector<LineSegment*>myEdges;
	vector<LineSegment*>yourEdges;

	auto l1 = f->findEdgePointer(t1,t2);
	auto l2 = f->findEdgePointer(t2,t3);
	auto l3 = f->findEdgePointer(t1,t3);

	myEdges.push_back(l1);
	myEdges.push_back(l2);
	myEdges.push_back(l3);

	auto s1 = (f->findEdgePointer(ls->t1,ls->t2));
	auto s2 = (f->findEdgePointer(ls->t2,ls->t3));
	auto s3 = (f->findEdgePointer(ls->t1,ls->t3));

	yourEdges.push_back(s1);
	yourEdges.push_back(s2);
	yourEdges.push_back(s3);

	bothEdges.insert(l1);
	bothEdges.insert(l2);
	bothEdges.insert(l3);
	bothEdges.insert(s1);
	bothEdges.insert(s2);
	bothEdges.insert(s3);

	if (bothEdges.size() != 5)
	{
		cerr<<"The triangle T and I do not share an edge.!"<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	else
	{
		for (auto itMy = myEdges.begin(); itMy != myEdges.end(); ++ itMy)
		{
			LineSegment* meE = *itMy;
			auto itf = find(yourEdges.begin(),yourEdges.end(),meE);
			if (itf != yourEdges.end())
			{
				sharingEdge = meE;
				break;
			}
		}
	}

	return sharingEdge;

}
Triangle* Triangle::findAItypeNeighborInFrontier(const int itype, vector<Triangle*> &frontier)
{
	Triangle* frontierNeighbor = NULL;
	if (frontier.empty())
	{
		cerr<<"The frontier should not be empty"<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	if (iTypeTriangleNeighbors[itype].empty())
	{
		cerr<<"The TriangleNeighbors should not be empty"<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	for (auto itTri = iTypeTriangleNeighbors[itype].begin(); itTri != iTypeTriangleNeighbors[itype].end(); ++ itTri)
	{
		Triangle* Tri = *itTri;

		auto it = find(frontier.begin(),frontier.end(), Tri);

		if (it != frontier.end())
		{
			frontierNeighbor = *it;
			break;
		}

	}

	return frontierNeighbor;

}
Triangle* Triangle::findANeighborInFrontier(vector<Triangle*> &frontier)
{
	Triangle* frontierNeighbor = NULL;
	if (frontier.empty())
	{
		cerr<<"The frontier should not be empty"<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	if (TriangleNeighbors.empty())
	{
		cerr<<"The TriangleNeighbors should not be empty"<<endl;
		cerr<<"Error in "<<__FILE__<<" at "<<__LINE__<<" in "<<__FUNCTION__<<endl;
		exit(EXIT_FAILURE);
	}
	for (auto itTri = TriangleNeighbors.begin(); itTri != TriangleNeighbors.end(); ++ itTri)
	{
		Triangle* Tri = *itTri;

		auto it = find(frontier.begin(),frontier.end(), Tri);

		if (it != frontier.end())
		{
			frontierNeighbor = *it;
			break;
		}

	}

	return frontierNeighbor;

}
Triangle::Triangle()
{
	Tri_ID = 0;
	Tri_field = NULL;
	t1 = NULL;
	t2 = NULL;
	t3 = NULL;
	Normalx = 0;
	Normaly = 0;
	Normalz = 0;
	triImportance = 0;
	broadcastCovered = 0;

	v2vMap.clear();
	v2theta.clear();
}

Triangle::~Triangle()
{
}

bool Triangle::ContainSensor(Sensor *s)
{
	int count = 0;
	if(s->CanFindASensorLandmarkNeighbors(t1))
	{
		++count;
	}
	if(s->CanFindASensorLandmarkNeighbors(t2))
	{
		++count;
	}
	if(s->CanFindASensorLandmarkNeighbors(t3))
	{
		++count;
	}
	if(count >= 3)
		return true;
	else
		return false;
}

void Triangle::iTypeEstablishTriNeighborhood(const int itype)
{
	for(unsigned int i = 0; i < Tri_field->iTypeTriangleSet[itype].size(); ++i)
	{
		Triangle *t = &(Tri_field->iTypeTriangleSet[itype][i]);
		if(this != t)
		{
			if(    (this->t1 == t->t1 && this->t2 == t->t2 && this->t3 != t->t3) 
				|| (this->t1 == t->t2 && this->t2 == t->t1 && this->t3 != t->t3) 
				|| (this->t1 == t->t1 && this->t2 == t->t3 && this->t3 != t->t2) 
				|| (this->t1 == t->t3 && this->t2 == t->t1 && this->t3 != t->t2) 
				|| (this->t1 == t->t3 && this->t2 == t->t2 && this->t3 != t->t1)
				|| (this->t1 == t->t2 && this->t2 == t->t3 && this->t3 != t->t1)

				|| (this->t3 == t->t1 && this->t1 == t->t2 && this->t2 != t->t3) 
				|| (this->t3 == t->t2 && this->t1 == t->t1 && this->t2 != t->t3) 
				|| (this->t3 == t->t1 && this->t1 == t->t3 && this->t2 != t->t2) 
				|| (this->t3 == t->t3 && this->t1 == t->t1 && this->t2 != t->t2) 
				|| (this->t3 == t->t3 && this->t1 == t->t2 && this->t2 != t->t1)
				|| (this->t3 == t->t2 && this->t1 == t->t3 && this->t2 != t->t1)

				|| (this->t2 == t->t1 && this->t3 == t->t2 && this->t1 != t->t3) 
				|| (this->t2 == t->t2 && this->t3 == t->t1 && this->t1 != t->t3) 
				|| (this->t2 == t->t1 && this->t3 == t->t3 && this->t1 != t->t2) 
				|| (this->t2 == t->t3 && this->t3 == t->t1 && this->t1 != t->t2) 
				|| (this->t2 == t->t3 && this->t3 == t->t2 && this->t1 != t->t1)
				|| (this->t2 == t->t2 && this->t3 == t->t3 && this->t1 != t->t1))
			{
				this->iTypeTriangleNeighbors[itype].insert(t);
			}
		}
	}
}
void Triangle::EstablishNeighborhood()
{
	FILE *fp;
	errno_t err;
	err = fopen_s(&fp, "triangle_vicility", "a+");
	if(!err)
	{
		fprintf(fp, "%d ", this->Tri_ID);		//Triangle ID;
		for(unsigned int i = 0; i < Tri_field->TriangleSet.size(); ++i)
		{
			Triangle *t = &Tri_field->TriangleSet[i];
			if(this != t)
			{
				if( (this->t1 == t->t1 && this->t2 == t->t2 && this->t3 != t->t3) 
					|| (this->t1 == t->t2 && this->t2 == t->t1 && this->t3 != t->t3) 
					|| (this->t1 == t->t1 && this->t2 == t->t3 && this->t3 != t->t2) 
					|| (this->t1 == t->t3 && this->t2 == t->t1 && this->t3 != t->t2) 
					|| (this->t1 == t->t3 && this->t2 == t->t2 && this->t3 != t->t1)
					|| (this->t1 == t->t2 && this->t2 == t->t3 && this->t3 != t->t1)
					
					|| (this->t3 == t->t1 && this->t1 == t->t2 && this->t2 != t->t3) 
					|| (this->t3 == t->t2 && this->t1 == t->t1 && this->t2 != t->t3) 
					|| (this->t3 == t->t1 && this->t1 == t->t3 && this->t2 != t->t2) 
					|| (this->t3 == t->t3 && this->t1 == t->t1 && this->t2 != t->t2) 
					|| (this->t3 == t->t3 && this->t1 == t->t2 && this->t2 != t->t1)
					|| (this->t3 == t->t2 && this->t1 == t->t3 && this->t2 != t->t1)
					
					|| (this->t2 == t->t1 && this->t3 == t->t2 && this->t1 != t->t3) 
					|| (this->t2 == t->t2 && this->t3 == t->t1 && this->t1 != t->t3) 
					|| (this->t2 == t->t1 && this->t3 == t->t3 && this->t1 != t->t2) 
					|| (this->t2 == t->t3 && this->t3 == t->t1 && this->t1 != t->t2) 
					|| (this->t2 == t->t3 && this->t3 == t->t2 && this->t1 != t->t1)
					|| (this->t2 == t->t2 && this->t3 == t->t3 && this->t1 != t->t1))
				{
					this->TriangleNeighbors.insert(t);

					fprintf(fp, "%d ", t->Tri_ID);		//Triangle neighbor ID;
					
				}
			}
		}
		fprintf(fp, "\n");
	}
	else
	{
		printf("Can not open file 'triangle_vicility'\n");
		exit(0);
	}
	fclose(fp);
}

#if 0
void Triangle::EstablishNeighborhood()
{
	for(unsigned int i = 0; i < Tri_field->TriangleSet.size(); ++i)
	{
		Triangle *t = &Tri_field->TriangleSet[i];
		if(this != t)
		{
			if( (this->t1 == t->t1 && this->t2 == t->t2 && this->t3 != t->t3) 
				|| (this->t1 == t->t2 && this->t2 == t->t1 && this->t3 != t->t3) 
				|| (this->t1 == t->t1 && this->t2 == t->t3 && this->t3 != t->t2) 
				|| (this->t1 == t->t3 && this->t2 == t->t1 && this->t3 != t->t2) 
				|| (this->t1 == t->t3 && this->t2 == t->t2 && this->t3 != t->t1)
				|| (this->t1 == t->t2 && this->t2 == t->t3 && this->t3 != t->t1)
				|| (this->t3 == t->t1 && this->t1 == t->t2 && this->t2 != t->t3) 
				|| (this->t3 == t->t2 && this->t1 == t->t1 && this->t2 != t->t3) 
				|| (this->t3 == t->t1 && this->t1 == t->t3 && this->t2 != t->t2) 
				|| (this->t3 == t->t3 && this->t1 == t->t1 && this->t2 != t->t2) 
				|| (this->t3 == t->t3 && this->t1 == t->t2 && this->t2 != t->t1)
				|| (this->t3 == t->t2 && this->t1 == t->t3 && this->t2 != t->t1)
				|| (this->t2 == t->t1 && this->t3 == t->t2 && this->t1 != t->t3) 
				|| (this->t2 == t->t2 && this->t3 == t->t1 && this->t1 != t->t3) 
				|| (this->t2 == t->t1 && this->t3 == t->t3 && this->t1 != t->t2) 
				|| (this->t2 == t->t3 && this->t3 == t->t1 && this->t1 != t->t2) 
				|| (this->t2 == t->t3 && this->t3 == t->t2 && this->t1 != t->t1)
				|| (this->t2 == t->t2 && this->t3 == t->t3 && this->t1 != t->t1))
			{
				this->TriangleNeighbors.insert(t);
			}
		}
	}
}
#endif

void Triangle::GetNormalVector()
{
	//Get normal of the triangle
	double N1x = this->t1->location.x - this->t2->location.x;
	double N1y = this->t1->location.y - this->t2->location.y;
	double N1z = this->t1->location.z - this->t2->location.z;
	double N2x = this->t1->location.x - this->t3->location.x;
	double N2y = this->t1->location.y - this->t3->location.y;
	double N2z = this->t1->location.z - this->t3->location.z;
	//Normal Vector
	double Nx = N1y/N1x*(N1x*N2z - N2x*N1z)/(N1x*N2y - N2x*N1y) - N1z/N1x;
	double Ny = (N2x*N1z - N1x*N2z)/(N1x*N2y  - N2x*N1y);
	double Nz = 1;
	if (N1z == 0 || N1x == 0 || N1y == 0 || N2z == 0 || N2x == 0 || N2y == 0 || (N1x*N2y-N2x*N1y)==0)
	{
		fprintf_s(stdout, "Vector problem!\n");
		//exit(0);
		if(N1x == 0)
			N1x = 0.01;
		if(N1y == 0)
			N1y = 0.01;
		if(N1z == 0)
			N1z = 0.01;
		if(N2x == 0)
			N2x = 0.01;
		if(N2y == 0)
			N2y = 0.01;
		if(N2z == 0)
			N2z = 0.01;
		if((N1x*N2y-N2x*N1y)==0)
			N1x += 0.01;
	}
	double Modulo = sqrt(Nx*Nx + Ny*Ny + Nz*Nz);
	Normalx = Nx/Modulo;
	Normaly = Ny/Modulo;
	Normalz = Nz/Modulo;
}

//////////////////////////////////////////////////////////////////////////Hole//////////////////////////////////////////////////////////////////////////
Hole::Hole(int holeID_, std::set<Sensor*> HoleSensorSet_)
{
	holeID = holeID_;
	HoleSensorSet.swap(HoleSensorSet_);
}

Hole::Hole()
{
	holeID = 0xffff;
	HoleSensorSet.clear();
}

Hole::~Hole()
{
	HoleSensorSet.clear();
}
						