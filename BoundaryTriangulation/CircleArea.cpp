#include "geometry3D.h"
#include "sensor.h"
#include "common.h"
#include <algorithm>
#include "assert.h"
#include "random.h"

CircleArea::CircleArea()
{
	Center = NULL;
	Circle.clear();
	AreaSensorSet.clear();
	Radiu = 0;
	contain_convex_number = 0;
	AvgCriticality = 0;
}

CircleArea::CircleArea(Field *_field, Sensor *_center, int _radiu)
{
	Center = _center;
	Radiu = _radiu;
	cirfield = _field;
	contain_convex_number = 0;
	int depth = 0;
	int out_circle_number = 0;		//number of sensors which is not within_circle
	double sum_criticality = 0;
	set<Sensor*> *frontier = new set<Sensor*>;
	Center->broadcastCovered = 1;
	if(Center->convex)
	{
		++contain_convex_number;
	}
	Center->within_circle = true;
	Center->asCenter = true;
	frontier->insert( Center );
	AreaSensorSet.insert(Center);
	while(depth != Radiu && !frontier->empty())
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
					if(s->convex)
					{
						++contain_convex_number;
					}
				}
			}
		}
		if(!newFrontier->empty())
		{
			++depth;
		}

		//dealing newFrontier
		set<Sensor*>::iterator newFrontierpos;

		if(depth == Radiu)
		{
			for(newFrontierpos = newFrontier->begin(); newFrontierpos != newFrontier->end(); ++newFrontierpos)
			{
				Sensor *s = *newFrontierpos;
				s->broadcastCovered = 1;
				if(s->onEdge_Degree != 0)
				{
					if(!s->within_circle)
					{
						++out_circle_number;
						s->within_circle = true;
					}
					s->on_boundary_circle = true;
					Circle.insert( make_pair(s->Criticality, s) );
					sum_criticality += s->Criticality;
					AreaSensorSet.insert(s);
				}
			}
		}
		else
		{
			for(newFrontierpos = newFrontier->begin(); newFrontierpos != newFrontier->end(); ++newFrontierpos)
			{
				Sensor *s = *newFrontierpos;
				s->broadcastCovered = 1;
				if(s->onEdge_Degree != 0)
				{
					if(!s->within_circle)
					{
						++out_circle_number;
						s->within_circle = true;
					}
					sum_criticality += s->Criticality;
					AreaSensorSet.insert(s);
				}
			}
		}

		delete frontier;
		frontier = newFrontier;
	}
	delete frontier;
	AvgCriticality = sum_criticality/(double)AreaSensorSet.size();
	for(int i = 0; i < cirfield->nSensors; ++i)
	{
		cirfield->sensorPool[i].broadcastCovered = 0;
	}
	//if(out_circle_number < cirfield->avg1HopNeighborhoodSize)
	//{
	//	Circle.clear();
	//	printf("Clear a circle!\n");
	//}
}

CircleArea::~CircleArea()
{
}

multimap<double, Sensor*, less<double> > CircleArea::GetCircle()
{
	return Circle;
}

void CircleArea::ComputingCircleAvgHopDist()
{
	multimap<double, Sensor*, less<double> >::iterator Circlepos1;
	multimap<double, Sensor*, less<double> >::iterator Circlepos2;
	double sum_dist_circle = 0;
	for(Circlepos1 = Circle.begin(); Circlepos1 != Circle.end(); ++Circlepos1)
	{
		Sensor *start = Circlepos1->second;
		double sum_dist_sensor = 0;
		double avg_dist_sensor = 0;
		int bad_hop = 0;
		for(Circlepos2 = Circle.begin(); Circlepos2 != Circle.end(); ++Circlepos2)
		{
			Sensor *end = Circlepos2->second;
			if(start->index != end->index)
			{
				if(start->UnknowDistTo(end))
				{
					int hop = start->FloodingforGetHopsBetweenStartandEnd(end);
					if(hop<=2)
					{
						++bad_hop;
					}
					sum_dist_sensor += hop;
					start->KnownHopDistEndID.push_back(end->index);
					end->KnownHopDistEndID.push_back(start->index);
					start->ID_HopDist_List.insert(make_pair(end->index, hop));
					end->ID_HopDist_List.insert(make_pair(start->index, hop));
				}
				else
				{
					for(map<int, int>::iterator IDHLpos = start->ID_HopDist_List.begin(); IDHLpos != start->ID_HopDist_List.end(); ++IDHLpos)
					{
						int IDH = IDHLpos->first;
						if(end->index == IDH)
						{
							int h = IDHLpos->second;
							if(h<=2)
							{
								++bad_hop;
							}
							sum_dist_sensor += h;
						}
					}
				}
			}
		}
		avg_dist_sensor = sum_dist_sensor/(Circle.size()-1-bad_hop);
		sum_dist_circle += avg_dist_sensor;
	}
	AvgHopdist = sum_dist_circle/Circle.size();
	for(multimap<double, Sensor*, less<double> >::iterator Circlepos = Circle.begin(); Circlepos != Circle.end(); ++Circlepos)
	{
		Sensor *c = Circlepos->second;
		c->KnownHopDistEndID.clear();
		c->ID_HopDist_List.clear();
	}
}