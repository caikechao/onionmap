#include "geometry3D.h"
#include "sensor.h"
#include "common.h"
#include <algorithm>
#include "assert.h"
#include "random.h"

void Path::EstablishPath(Sensor *s1, Sensor *s2)
{
	start = s1;
	end = s2;
	path = start->PathFlood(end);
}

void Path::EstablishDirectPath(Sensor *s1, Sensor *s2)
{
	start = s1;
	end = s2;
	path = start->PathFloodforDirectPath(end);
}

void Path::EstablishPathBySideStep(Sensor *s1, Sensor *s2, set<Sensor*> area_sensor_set)
{
	start = s1;
	end = s2;
	path.clear();
	if(Path_Changed(start, area_sensor_set))
	{
		path.push_back(start);
	}
	else
	{
		printf("Fail to establish path by side-step!!!\n");
		//exit(0);	
	}
	for(int i = 0; i < pafield->nSensors; i++) 
	{
		pafield->sensorPool[i].broadcastCovered = 0;
		pafield->sensorPool[i].candidateParents.clear();
		pafield->sensorPool[i].parent = NULL;
		pafield->sensorPool[i].children.clear();
	}
}

void Path::EstablishLowestCriticalityShortestPath(Sensor *_start, Sensor *_end, set<Sensor*> area_sensor_set)
{
	start = _start;
	end = _end;
	path = start->PathCriticalityFlooding(end, area_sensor_set);
}

/*
void Path::FindConcavePoints(Sensor *_start, Sensor *_end)
{
	start = _start;
	end = _end;
	path = start->FindConcaveFlooding(end);
}*/


#if 0
bool Path::PathCrossing(Path *pa, int hopnum)
{
	int i;
	for(i = 0; i < path.size(); ++i)
	{
		if(path[i]->index == start->index || path[i]->index == end->index)
		{
			continue;
		}
		if(path[i]->NotNearEndsOfPath(this) && path[i]->SearchInKhopNeighbors(this, pa, hopnum))
		{
			return true;
		}
	}
	return false;
}
#endif

bool Path::Path_Changed(Sensor *floodsensor, set<Sensor*> area_sensor_set)
{
	floodsensor->broadcastCovered = 1;
	if(floodsensor == end)
	{
		return true;
	}
	else
	{
		bool newFrontier_empty = false;
		multimap<double, Sensor*, less<double> > *newFrontier = new multimap<double, Sensor*, less<double> >;
		for(SCITER neighborpos = floodsensor->neighbors.begin(); neighborpos != floodsensor->neighbors.end(); ++neighborpos) {		
			Sensor *n = *neighborpos;
			if (n->broadcastCovered == 0 /*&& n->onEdge_Degree > 0*/ && n->WithinArea(area_sensor_set)) {
				newFrontier->insert( make_pair(n->Criticality, n) );
				n->broadcastCovered = 1;
			}
		}
		//dealing with newFrontier
		if(newFrontier->empty())
		{
			newFrontier_empty = true;
		}
		for(multimap<double, Sensor*, less<double> >::iterator newFrontierpos = newFrontier->begin(); newFrontierpos != newFrontier->end(); ++newFrontierpos)
		{
			if(newFrontier->empty())
			{
				break;
			}
			else
			{
				floodsensor = newFrontierpos->second;
				if(Path_Changed(floodsensor, area_sensor_set))
				{
					path.push_back(floodsensor);
					delete newFrontier;
					return true;
				}
			}
		}
		delete newFrontier;
		return false;
	}
}

bool Path::PathCrossing(int hopnum)
{
	bool CancelTheLine = false;
	for(unsigned int p = 1; p < this->path.size()-1; ++p)
	{
		Sensor *lp = this->path[p];
		if(lp->BelongOtherIntersectLine(this->start->index, this->end->index))
		{
			CancelTheLine = true;
			//printf("Detect crossing!\n");
			break;
		}


		for(SCITER lpnpos = lp->neighbors.begin(); lpnpos != lp->neighbors.end(); ++lpnpos)
		{
			Sensor *lpn = *lpnpos;
			if(lpn->onEdge_Degree != 0 && lpn->BelongOtherIntersectLine(this->start->index, this->end->index))
			{
				CancelTheLine = true;
				//printf("Detect crossing!\n");
				break;
			}
		}
		if (CancelTheLine)
		{
			break;
		}
	}
	if(!CancelTheLine)
	{
		return false;
	}
	else
	{
		return true;
	}
}

bool Path::IntersectWith(LineSegment *dl, int hopnum)
{
	unsigned int i;
	unsigned int j;
	for(i = 0; i < path.size(); ++i)
	{
		if(path[i]->index == start->index || path[i]->index == end->index)
		{
			continue;
		}
		if(path[i]->NotNearEndsOfPath(this))
		{
			int depth = 0;
			set<Sensor*> *frontier = new set<Sensor*>;
			path[i]->broadcastCovered = 1;
			for(j = 0; j < dl->Lpath.path.size(); ++j)
			{
				Sensor *dlpathj = dl->Lpath.path[j];
				if(path[i]->index == dlpathj->index)
				{
					return true;		//intersect
				}
			}
			frontier->insert( path[i] );
			
			while(depth != hopnum && !frontier->empty())
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
					if(depth == hopnum)
					{
						for(j = 0; j < dl->Lpath.path.size(); ++j)
						{
							Sensor *dlpathj = dl->Lpath.path[j];
							if(s->index == dlpathj->index)
							{
								delete frontier;
								delete newFrontier;
								return true;		//intersect
							}
						}
					}
				}
				delete frontier;
				frontier = newFrontier;
			}
			delete frontier;
			for(j = 0; j < pafield->BoundarySensorSet.size(); j++) 
			{
				pafield->BoundarySensorSet[j]->broadcastCovered = 0;
			}
		}
		
	}
	return false;	//not intersect
}