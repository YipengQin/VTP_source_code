#ifndef GEODESIC_ALGORITHM_BASE
#define GEODESIC_ALGORITHM_BASE

#include "stdafx.h"
#include "geodesic_mesh.h"
#include "geodesic_algorithm_exact_elements.h"
#include "geodesic_constants_and_simple_functions.h"

namespace geodesic{

class GeodesicAlgorithmBase
{
public:

	GeodesicAlgorithmBase(geodesic::Mesh* mesh): m_mesh(mesh),
		m_edge_interval_lists_0(mesh->edges().size()),
		m_edge_interval_lists_1(mesh->edges().size())
	{
		// initialize statistics
		m_queue_max_size      = 0;
		m_windows_propagation = 0;
		m_windows_wavefront   = 0;
		m_windows_peak        = 0;

		// initialize window lists, similar to half-edge structure
		for (unsigned i = 0; i < m_edge_interval_lists_0.size(); ++i)
		{
			m_edge_interval_lists_0[i].initialize(&mesh->edges()[i]);
			m_edge_interval_lists_1[i].initialize(&mesh->edges()[i]);
		}
		
		std::queue<list_pointer> list_queue;
		edge_pointer edge = &(this->mesh()->edges()[0]);
		interval_list_0(edge)->start_vertex() = edge->v0();
		list_queue.push(interval_list_0(edge));

		while (!list_queue.empty())
		{
			list_pointer list = list_queue.front();
			list_queue.pop();

			edge_pointer edge = list->edge();
			face_pointer face;

			if (list == interval_list_0(edge))
				face = edge->adjacent_faces()[0];
			else
				face = edge->adjacent_faces().size() > 1 ? edge->adjacent_faces()[1] : NULL; // Boundary Case
			
			// opposite list
			if (list == interval_list_0(edge))
			{
				if (!interval_list_1(edge)->start_vertex())
				{
					interval_list_1(edge)->start_vertex() = edge->opposite_vertex(list->start_vertex());
					list_queue.push(interval_list_1(edge));
				}
			}
			else
			{
				if (!interval_list_0(edge)->start_vertex())
				{
					interval_list_0(edge)->start_vertex() = edge->opposite_vertex(list->start_vertex());
					list_queue.push(interval_list_0(edge));
				}
			}

			// neighbour list 1
			if (face)
			{
				edge_pointer edge_next;
				list_pointer list_next;

				edge_next = face->next_edge(edge, list->start_vertex());
				list_next = edge_next->adjacent_faces()[0] == face ? interval_list_0(edge_next) : interval_list_1(edge_next);

				if (!list_next->start_vertex())
				{
					list_next->start_vertex() = edge_next->opposite_vertex(list->start_vertex());
					list_queue.push(list_next);
				}

				// neighbour list 2
				edge_next = face->next_edge(edge, edge->opposite_vertex(list->start_vertex()));
				list_next = edge_next->adjacent_faces()[0] == face ? interval_list_0(edge_next) : interval_list_1(edge_next);

				if (!list_next->start_vertex())
				{
					list_next->start_vertex() = edge->opposite_vertex(list->start_vertex());
					list_queue.push(list_next);
				}
			}
		}

		// verify list links
		for (unsigned i = 0; i < this->mesh()->faces().size(); ++i)
		{
			face_pointer f = &(this->mesh()->faces()[i]);
			vertex_pointer v[3];
			for (unsigned j = 0; j < 3; ++j)
			{
				edge_pointer e = f->adjacent_edges()[j];

				if (e->adjacent_faces()[0] == f)
					v[j] = interval_list_0(e)->start_vertex();
				else
					v[j] = interval_list_1(e)->start_vertex();

				if ((!interval_list_0(e)->start_vertex()) || (!interval_list_1(e)->start_vertex()))
				{
					std::cout << "list link error" << std::endl;
					exit(1);
				}

				if (interval_list_0(e)->start_vertex() == interval_list_1(e)->start_vertex())
				{
					std::cout << "list link error" << std::endl;
					exit(1);
				}

				if (!((e->belongs(interval_list_0(e)->start_vertex())) && (e->belongs(interval_list_1(e)->start_vertex()))))
				{
					std::cout << "list link error" << std::endl;
					exit(1);
				}
			}
			if ((v[0] == v[1]) || (v[0] == v[2]) || (v[1] == v[2]))
			{
				std::cout << "list link error" << std::endl;
				exit(1);
			}
		}
	};

	virtual ~GeodesicAlgorithmBase(){};

	virtual void print_statistics()		//print info about timing and memory usage in the propagation step of the algorithm
	{
		std::cout << "propagation step took " << m_time_consumed << " seconds " << std::endl;
	};	

	geodesic::Mesh* mesh(){return m_mesh;};

	// propagate a window
	bool compute_propagated_parameters(double pseudo_x,
		double pseudo_y,
		double start,
		double end,		//start/end of the interval
		double alpha,	//corner angle
		double L,		//length of the new edge
		interval_pointer candidates,
		double d);		//if it is the last interval on the edge

	// intersection point on an edge
	double compute_positive_intersection(double start,
		double pseudo_x,
		double pseudo_y,
		double sin_alpha,
		double cos_alpha);

	inline bool calculate_triangle_parameters(list_pointer &list, Triangle &Tri); // calculate the parameters of the triangle to be propagated

	list_pointer interval_list_0(edge_pointer e)
	{
		return &m_edge_interval_lists_0[e->id()];
	};

	list_pointer interval_list_1(edge_pointer e)
	{
		return &m_edge_interval_lists_1[e->id()];
	};

protected:

	geodesic::Mesh* m_mesh;

	Triangle Tri; // the triangle to be propagated

	IntervalList wl_left, wl_right;

	std::vector<IntervalList> m_edge_interval_lists_0;		// windows propagated from adjacent_face[0] of the edge
	std::vector<IntervalList> m_edge_interval_lists_1;		// windows propagated from adjacent_face[1] of the edge

};

inline double GeodesicAlgorithmBase::compute_positive_intersection(double start,
	double pseudo_x,
	double pseudo_y,
	double sin_alpha,
	double cos_alpha)
{
	//assert(pseudo_y < 0);
	assert(pseudo_y <= 0);

	double denominator = sin_alpha*(pseudo_x - start) - cos_alpha*pseudo_y;
	if (denominator < 0.0)
	{
		return -1.0;
	}

	double numerator = -pseudo_y*start;

	if (numerator < 1e-30)
	{
		return 0.0;
	}

	if (denominator < 1e-30)
	{
		return -1.0;
	}

	return numerator / denominator;
}

inline bool GeodesicAlgorithmBase::compute_propagated_parameters(double pseudo_x,
	double pseudo_y,
	double begin,
	double end,		//start/end of the interval
	double alpha,	//corner angle
	double L,		//length of the new edge
	interval_pointer candidates,
	double d)
{
	assert(pseudo_y <= 0.0);
	assert(begin <= end);
	assert(begin >= 0);

	++m_windows_propagation; // Statistics

	interval_pointer p = candidates;

	double sin_alpha = sin(alpha);
	double cos_alpha = cos(alpha);

	//important: for the first_interval, this function returns zero only if the new edge is "visible" from the source
	//if the new edge can be covered only after turn_over, the value is negative (-1.0)
	double L1 = compute_positive_intersection(begin,
		pseudo_x,
		pseudo_y,
		sin_alpha,
		cos_alpha);

	if (L1 < 0 || L1 >= L) // Does not produce a window on the edge
		return false;

	double L2 = compute_positive_intersection(end,
		pseudo_x,
		pseudo_y,
		sin_alpha,
		cos_alpha);

	if (L2 < 0 || L2 >= L) // Covers vertex
	{
		p->start() = L1;
		p->stop() = L;
		p->pseudo_x() = cos_alpha*pseudo_x + sin_alpha*pseudo_y;
		p->pseudo_y() = -sin_alpha*pseudo_x + cos_alpha*pseudo_y;
		assert(p->pseudo_y() <= 0.0);

		return true;
	}
	else
	{
		// Does not cover vertex
		p->start() = L1;
		p->stop() = L2;
		p->pseudo_x() = cos_alpha*pseudo_x + sin_alpha*pseudo_y;
		p->pseudo_y() = -sin_alpha*pseudo_x + cos_alpha*pseudo_y;
		assert(p->pseudo_y() <= 0.0);

		return true;
	}
}

inline bool GeodesicAlgorithmBase::calculate_triangle_parameters(list_pointer &list, Triangle &Tri) // Calculate the parameters of the triangle to be propagated
{
	if (list->edge()->adjacent_faces().size() > 1)
	{
		Tri.bottom_edge = list->edge();

		if (list == interval_list_0(Tri.bottom_edge))
			Tri.face = Tri.bottom_edge->adjacent_faces()[1];
		else
			Tri.face = Tri.bottom_edge->adjacent_faces()[0];

		Tri.top_vertex = Tri.face->opposite_vertex(Tri.bottom_edge);
		Tri.left_vertex = list->start_vertex();
		Tri.right_vertex = Tri.bottom_edge->opposite_vertex(Tri.left_vertex);

		Tri.left_edge = Tri.face->next_edge(Tri.bottom_edge, Tri.left_vertex);
		Tri.right_edge = Tri.face->next_edge(Tri.bottom_edge, Tri.right_vertex);

		Tri.top_alpha = Tri.face->vertex_angle(Tri.top_vertex);
		Tri.left_alpha = Tri.face->vertex_angle(Tri.left_vertex);
		Tri.right_alpha = Tri.face->vertex_angle(Tri.right_vertex);

		if (Tri.left_edge->adjacent_faces()[0] == Tri.face)
			Tri.left_list = interval_list_0(Tri.left_edge);
		else
			Tri.left_list = interval_list_1(Tri.left_edge);

		if (Tri.right_edge->adjacent_faces()[0] == Tri.face)
			Tri.right_list = interval_list_0(Tri.right_edge);
		else
			Tri.right_list = interval_list_1(Tri.right_edge);

		return false;
	}
	else
	{
		return true;
	}
}


}//geodesic

#endif
