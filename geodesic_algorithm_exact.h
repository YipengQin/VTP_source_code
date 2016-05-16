#ifndef GEODESIC_ALGORITHM_EXACT
#define GEODESIC_ALGORITHM_EXACT

#include "stdafx.h"
#include "geodesic_memory.h"
#include "geodesic_algorithm_base.h"
#include "geodesic_algorithm_exact_elements.h"

namespace geodesic {

	class GeodesicAlgorithmExact : public GeodesicAlgorithmBase
	{
	public:

		// basic functions related to class
		GeodesicAlgorithmExact(geodesic::Mesh* mesh) :
			GeodesicAlgorithmBase(mesh),
			m_memory_allocator(mesh->edges().size(), mesh->edges().size())
		{};
		~GeodesicAlgorithmExact() {};
		void clear() { m_memory_allocator.clear(); };

		// main entry
		void propagate(unsigned source);

		// print the resulting statistics
		void print_statistics();

	private:

		// simple functions
		void initialize_propagation_data();
		void create_pseudo_source_windows(vertex_pointer &v, bool UpdateFIFOQueue);
		void erase_from_queue(vertex_pointer v);

		// propagate a windows list (Rule 1)
		void find_separating_point(list_pointer &list); // find the separating point of the windows and the list
		void propagate_windows_to_two_edges(list_pointer &list); // propagates windows to two edges accross a triangle face

		// pairwise windows checking (Rule 2)
		void check_with_vertices(list_pointer &list);
		windows_state check_between_two_windows(interval_pointer &w1, interval_pointer &w2); // Check two neighbouring crossing windows on same edge
		void pairwise_windows_checking(list_pointer &list); // Check crossing windows on same edge

		// main operation
		void propagate_one_windows_list(list_pointer &list);

		// member variables
		//std::multiset<vertex_pointer, Vertex> m_vertex_queue; // priority queue for vertices
		std::set<vertex_pointer, Vertex> m_vertex_queue;
		std::queue<list_pointer> m_list_queue;                // FIFO queue for lists
		MemoryAllocator<Interval> m_memory_allocator;		  // quickly allocate and deallocate intervals 

		unsigned m_source;
	};

	//----------------- simple functions ---------------------
	inline void GeodesicAlgorithmExact::initialize_propagation_data()
	{
		clear();

		// initialize source's parameters
		vertex_pointer source = &(this->mesh()->vertices()[m_source]);
		source->geodesic_distance() = 0;
		source->state() = Vertex::INSIDE;

		// initialize windows around source
		create_pseudo_source_windows(source, false);
	}

	inline void GeodesicAlgorithmExact::erase_from_queue(vertex_pointer v)
	{
		assert(m_vertex_queue.count(v) <= 1);

		std::multiset<vertex_pointer, Vertex>::iterator it = m_vertex_queue.find(v);
		if (it != m_vertex_queue.end())
			m_vertex_queue.erase(it);
	}

	inline void GeodesicAlgorithmExact::create_pseudo_source_windows(vertex_pointer &pseudo_source, bool inside_traversed_area)
	{
		// update vertices around pseudo_source
		for (unsigned i = 0; i < pseudo_source->adjacent_edges().size(); ++i)
		{
			edge_pointer   edge_it = pseudo_source->adjacent_edges()[i];
			vertex_pointer vert_it = edge_it->opposite_vertex(pseudo_source);

			double distance = pseudo_source->geodesic_distance() + edge_it->length();
			if (distance < vert_it->geodesic_distance())
			{
				m_vertex_queue.erase(vert_it);

				vert_it->geodesic_distance() = distance;
				if (vert_it->state() == Vertex::OUTSIDE)
					vert_it->state() = Vertex::FRONT;

				vert_it->incident_face() = edge_it->adjacent_faces()[0];
				edge_pointer next_edge = vert_it->incident_face()->next_edge(edge_it, pseudo_source);
				vert_it->incident_point() = (next_edge->v0() == pseudo_source) ? 0 : next_edge->length();

				m_vertex_queue.insert(vert_it);
			}
		}

		// update pseudo_source windows around pseudo_source
		for (unsigned i = 0; i < pseudo_source->adjacent_faces().size(); ++i)
		{
			face_pointer face_it = pseudo_source->adjacent_faces()[i];
			edge_pointer edge_it = face_it->opposite_edge(pseudo_source);
			list_pointer list = (edge_it->adjacent_faces()[0] == face_it) ? list = interval_list_0(edge_it) : list = interval_list_1(edge_it);

			// create a window
			interval_pointer candidate = new Interval;

			candidate->start() = 0;
			candidate->stop() = edge_it->length();
			candidate->d() = pseudo_source->geodesic_distance();
			double angle = face_it->vertex_angle(list->start_vertex());
			double length = face_it->next_edge(edge_it, list->start_vertex())->length();
			candidate->pseudo_x() = cos(angle) * length;
			candidate->pseudo_y() = -sin(angle) * length;

			// insert into list
			list->push_back(candidate);

			// push into M_LIST_QUEUE if inside traversed area
			if ((inside_traversed_area) && ((edge_it->v0()->state() != Vertex::FRONT) || (edge_it->v1()->state() != Vertex::FRONT)))
				m_list_queue.push(list);

			// Statistics
			++m_windows_wavefront;
			if (m_windows_peak < m_windows_wavefront)
				m_windows_peak = m_windows_wavefront;

		}
	}

	//----------------- propagate a windows list (Rule 1) ---------------------
	inline void GeodesicAlgorithmExact::find_separating_point(list_pointer &list)
	{
		const double LOCAL_EPSILON = 1e-20 * list->edge()->length(); // numerical issue

		double L = Tri.left_edge->length();
		double top_x = L * cos(Tri.left_alpha);
		double top_y = L * sin(Tri.left_alpha);

		Vertex top_t; // temporary top_vertex
		memcpy(&top_t, Tri.top_vertex, sizeof(Vertex));
		top_t.geodesic_distance() = GEODESIC_INF;

		interval_pointer iter = list->begin();

		double wlist_sp = 0;
		double wlist_pseudo_x = 0;
		double wlist_pseudo_y = 0;

		while (iter != NULL)
		{
			interval_pointer &w = iter;

			double w_sp = w->pseudo_x() - w->pseudo_y() * ((top_x - w->pseudo_x()) / (top_y - w->pseudo_y()));
			double distance = GEODESIC_INF;

			// shortest path from the window
			if ((w_sp - w->start() > LOCAL_EPSILON) && (w_sp - w->stop() < -LOCAL_EPSILON))
			{
				distance = w->d() + sqrt((top_x - w->pseudo_x()) * (top_x - w->pseudo_x()) + (top_y - w->pseudo_y()) * (top_y - w->pseudo_y()));
				w->shortest_distance() = distance;
			}
			else if (w_sp - w->start() <= LOCAL_EPSILON)
			{
				distance = w->d() + sqrt((top_x - w->start()) * (top_x - w->start()) + top_y * top_y) + sqrt((w->start() - w->pseudo_x()) * (w->start() - w->pseudo_x()) + w->pseudo_y() * w->pseudo_y());
				w->shortest_distance() = distance;
				w_sp = w->start();
			}
			else if (w_sp - w->stop() >= -LOCAL_EPSILON)
			{
				distance = w->d() + sqrt((top_x - w->stop()) * (top_x - w->stop()) + top_y * top_y) + sqrt((w->stop() - w->pseudo_x()) * (w->stop() - w->pseudo_x()) + w->pseudo_y() * w->pseudo_y());
				w->shortest_distance() = distance;
				w_sp = w->stop();
			}

			// update information at top_t
			if (distance < top_t.geodesic_distance())
			{
				top_t.geodesic_distance() = distance;
				top_t.incident_face() = Tri.face;
				top_t.incident_point() = (list->start_vertex() == list->edge()->v0()) ? w_sp : list->edge()->length() - w_sp;
				wlist_sp = w_sp;
				wlist_pseudo_x = w->pseudo_x();
				wlist_pseudo_y = w->pseudo_y();
			}
			w->sp() = w_sp;

			iter = iter->next();
		}

		// update top_vertex and M_VERTEX_QUEUE
		if (top_t.geodesic_distance() < Tri.top_vertex->geodesic_distance())
		{
			if (Tri.top_vertex->state() == Vertex::FRONT) erase_from_queue(Tri.top_vertex);
			memcpy(Tri.top_vertex, &top_t, sizeof(Vertex));
			if (Tri.top_vertex->state() == Vertex::FRONT) m_vertex_queue.insert(Tri.top_vertex);

			if ((Tri.top_vertex->state() == Vertex::INSIDE) && (Tri.top_vertex->saddle_or_boundary()))
				create_pseudo_source_windows(Tri.top_vertex, true); // handle saddle vertex
		}

		list->sp() = wlist_sp;
		list->pseudo_x() = wlist_pseudo_x;
		list->pseudo_y() = wlist_pseudo_y;
	}

	inline void GeodesicAlgorithmExact::propagate_windows_to_two_edges(list_pointer &list)
	{
		const double LOCAL_EPSILON = 1e-8 * list->edge()->length(); // numerical issue

		interval_pointer iter = list->begin();
		interval_pointer iter_t;

		enum PropagationDirection
		{
			LEFT,
			RIGHT,
			BOTH
		};

		PropagationDirection direction;

		while (!list->empty() && (iter != NULL))
		{
			interval_pointer &w = iter;
			assert(w->start() <= w->stop());

			if (w->sp() < list->sp() - LOCAL_EPSILON)
			{
				// only propagate to left edge
				double Intersect_X, Intersect_Y;

				// judge the positions of the two windows
				CalculateIntersectionPoint(list->pseudo_x(), list->pseudo_y(), list->sp(), 0, w->pseudo_x(), w->pseudo_y(), w->stop(), 0, Intersect_X, Intersect_Y);
				if ((w->stop() < list->sp()) || ((Intersect_Y <= 0) && (Intersect_Y >= list->pseudo_y()) && (Intersect_Y >= w->pseudo_y())))
				{
					direction = PropagationDirection::LEFT;
				}
				else
				{
					direction = PropagationDirection::BOTH;
				}				
			}
			else if (w->sp() > list->sp() + LOCAL_EPSILON)
			{
				// only propagate to right edge
				double Intersect_X, Intersect_Y;

				// judge the positions of the two windows
				CalculateIntersectionPoint(list->pseudo_x(), list->pseudo_y(), list->sp(), 0, w->pseudo_x(), w->pseudo_y(), w->start(), 0, Intersect_X, Intersect_Y);
				if ((w->start() > list->sp())||((Intersect_Y <= 0) && (Intersect_Y >= list->pseudo_y()) && (Intersect_Y >= w->pseudo_y())))
				{
					direction = PropagationDirection::RIGHT;
				}
				else
				{
					direction = PropagationDirection::BOTH;
				}	
			}
			else
			{
				// propagate to both edges
				direction = PropagationDirection::BOTH;
			}

			bool ValidPropagation;
			interval_pointer right_w;

			switch (direction) {
			case PropagationDirection::LEFT:
				ValidPropagation = compute_propagated_parameters(w->pseudo_x(),
					w->pseudo_y(),
					w->start(),
					w->stop(),
					Tri.left_alpha,
					Tri.left_edge->length(),
					w,
					w->d());

				iter_t = iter->next();
				if (ValidPropagation)
				{
					list->erase(w);
					wl_left.push_back(w);
				}
				else
				{
					list->erase(w);
					delete w;
					--m_windows_wavefront;
				}
				iter = iter_t;
				break;

			case PropagationDirection::RIGHT:
				ValidPropagation = compute_propagated_parameters(Tri.bottom_edge->length() - w->pseudo_x(),
					w->pseudo_y(),
					Tri.bottom_edge->length() - w->stop(),
					Tri.bottom_edge->length() - w->start(),
					Tri.right_alpha,
					Tri.right_edge->length(),
					w,
					w->d());

				iter_t = iter->next();
				if (ValidPropagation)
				{
					double length = Tri.right_edge->length(); // invert window
					double start = length - w->stop();
					w->stop() = length - w->start();
					w->start() = start;
					w->pseudo_x() = length - w->pseudo_x();

					list->erase(w);
					wl_right.push_back(w);
				}
				else
				{
					list->erase(w);
					delete w;
					--m_windows_wavefront;
				}
				iter = iter_t;
				break;

			case PropagationDirection:: BOTH:
				right_w = new Interval;
				memcpy(right_w, w, sizeof(Interval));

				ValidPropagation = compute_propagated_parameters(w->pseudo_x(),
					w->pseudo_y(),
					w->start(),
					w->stop(),
					Tri.face->vertex_angle(Tri.left_vertex),
					Tri.left_edge->length(),
					w,
					w->d());

				iter_t = iter->next();
				if (ValidPropagation)
				{
					list->erase(w);
					wl_left.push_back(w);
				}
				else
				{
					list->erase(w);
					delete w;
					--m_windows_wavefront;
				}
				iter = iter_t;

				ValidPropagation = compute_propagated_parameters(Tri.bottom_edge->length() - right_w->pseudo_x(),
					right_w->pseudo_y(),
					Tri.bottom_edge->length() - right_w->stop(),
					Tri.bottom_edge->length() - right_w->start(),
					Tri.face->vertex_angle(Tri.right_vertex),
					Tri.right_edge->length(),
					right_w,
					right_w->d());

				if (ValidPropagation)
				{
					// invert window
					double length = Tri.right_edge->length();
					double start = length - right_w->stop();
					right_w->stop() = length - right_w->start();
					right_w->start() = start;
					right_w->pseudo_x() = length - right_w->pseudo_x();

					wl_right.push_back(right_w);

					++m_windows_wavefront;
					if (m_windows_peak < m_windows_wavefront)
						m_windows_peak = m_windows_wavefront;
				}
				else
				{
					delete right_w;
				}
				break;

			default:
				break;
			}
		}
	}

	//----------------- pairwise windows checking (Rule 2) ----------------------
	inline void GeodesicAlgorithmExact::check_with_vertices(list_pointer &list)
	{
		if (list->empty()) return;

		interval_pointer iter = list->begin();
		interval_pointer iter_t;

		while ((!list->empty()) && (iter != NULL))
		{
			interval_pointer &w = iter;
			bool w_survive = true;

			edge_pointer   e = list->edge();
			vertex_pointer v1 = list->start_vertex();
			vertex_pointer v2 = e->opposite_vertex(v1);
			double d1 = GEODESIC_INF;

			d1 = w->d() + sqrt((w->stop() - w->pseudo_x()) * (w->stop() - w->pseudo_x()) + w->pseudo_y() * w->pseudo_y());
			if (v1->geodesic_distance() + w->stop() < d1)
				w_survive = false;

			d1 = w->d() + sqrt((w->start() - w->pseudo_x()) * (w->start() - w->pseudo_x()) + w->pseudo_y() * w->pseudo_y());
			if (v2->geodesic_distance() + e->length() - w->start() < d1)
				w_survive = false;


			iter_t = iter;
			iter = iter->next();

			if (!w_survive)
			{
				list->erase(iter_t);
				delete iter_t;
				--m_windows_wavefront;
			}
		}
	}

	inline windows_state GeodesicAlgorithmExact::check_between_two_windows(interval_pointer &w1, interval_pointer &w2)
	{
		double NUMERCIAL_EPSILON = 1 - 1e-12;
		// we implement the discussed 6 cases as follows for simplicity

		if ((w1->start() >= w2->start()) && (w1->start() <= w2->stop())) // w1->start
		{
			double Intersect_X, Intersect_Y;

			// judge the order of the two windows
			CalculateIntersectionPoint(w2->pseudo_x(), w2->pseudo_y(), w1->start(), 0, w1->pseudo_x(), w1->pseudo_y(), w1->stop(), 0, Intersect_X, Intersect_Y);

			if ((Intersect_Y <= 0) && (Intersect_Y >= w1->pseudo_y()) && (Intersect_Y >= w2->pseudo_y()))
			{
				double d1, d2;
				d1 = w1->d() + sqrt((w1->start() - w1->pseudo_x()) * (w1->start() - w1->pseudo_x()) + (w1->pseudo_y()) * (w1->pseudo_y()));
				d2 = w2->d() + sqrt((w1->start() - w2->pseudo_x()) * (w1->start() - w2->pseudo_x()) + (w2->pseudo_y()) * (w2->pseudo_y()));

				if (d2 < d1 * NUMERCIAL_EPSILON)
					return w1_invalid;
				if (d1 < d2 * NUMERCIAL_EPSILON)
					w2->start() = w1->start();
			}
		}

		if ((w1->stop() >= w2->start()) && (w1->stop() <= w2->stop())) // w1->stop
		{
			double Intersect_X, Intersect_Y;

			// judge the order of the two windows
			CalculateIntersectionPoint(w2->pseudo_x(), w2->pseudo_y(), w1->stop(), 0, w1->pseudo_x(), w1->pseudo_y(), w1->start(), 0, Intersect_X, Intersect_Y);

			if ((Intersect_Y <= 0) && (Intersect_Y >= w1->pseudo_y()) && (Intersect_Y >= w2->pseudo_y()))
			{
				double d1, d2;
				d1 = w1->d() + sqrt((w1->stop() - w1->pseudo_x()) * (w1->stop() - w1->pseudo_x()) + (w1->pseudo_y()) * (w1->pseudo_y()));
				d2 = w2->d() + sqrt((w1->stop() - w2->pseudo_x()) * (w1->stop() - w2->pseudo_x()) + (w2->pseudo_y()) * (w2->pseudo_y()));

				if (d2 < d1 * NUMERCIAL_EPSILON)
					return w1_invalid;
				if (d1 < d2 * NUMERCIAL_EPSILON)
					w2->stop() = w1->stop();
			}
		}

		if ((w2->start() >= w1->start()) && (w2->start() <= w1->stop())) // w2->start
		{
			double Intersect_X, Intersect_Y;

			// judge the previous order of the two windows
			CalculateIntersectionPoint(w1->pseudo_x(), w1->pseudo_y(), w2->start(), 0, w2->pseudo_x(), w2->pseudo_y(), w2->stop(), 0, Intersect_X, Intersect_Y);

			if ((Intersect_Y <= 0) && (Intersect_Y >= w1->pseudo_y()) && (Intersect_Y >= w2->pseudo_y()))
			{
				double d1, d2;
				d1 = w1->d() + sqrt((w2->start() - w1->pseudo_x()) * (w2->start() - w1->pseudo_x()) + (w1->pseudo_y()) * (w1->pseudo_y()));
				d2 = w2->d() + sqrt((w2->start() - w2->pseudo_x()) * (w2->start() - w2->pseudo_x()) + (w2->pseudo_y()) * (w2->pseudo_y()));

				if (d1 < d2 * NUMERCIAL_EPSILON)
					return w2_invalid;
				if (d2 < d1 * NUMERCIAL_EPSILON)
					w1->start() = w2->start();
			}
		}

		if ((w2->stop() >= w1->start()) && (w2->stop() <= w1->stop())) // w2->stop
		{
			double Intersect_X, Intersect_Y;

			// judge the previous order of the two windows
			CalculateIntersectionPoint(w1->pseudo_x(), w1->pseudo_y(), w2->stop(), 0, w2->pseudo_x(), w2->pseudo_y(), w2->start(), 0, Intersect_X, Intersect_Y);

			if ((Intersect_Y <= 0) && (Intersect_Y >= w1->pseudo_y()) && (Intersect_Y >= w2->pseudo_y()))
			{
				double d1, d2;
				d1 = w1->d() + sqrt((w2->stop() - w1->pseudo_x()) * (w2->stop() - w1->pseudo_x()) + (w1->pseudo_y()) * (w1->pseudo_y()));
				d2 = w2->d() + sqrt((w2->stop() - w2->pseudo_x()) * (w2->stop() - w2->pseudo_x()) + (w2->pseudo_y()) * (w2->pseudo_y()));

				if (d1 < d2 * NUMERCIAL_EPSILON)
					return w2_invalid;
				if (d2 < d1 * NUMERCIAL_EPSILON)
					w1->stop() = w2->stop();
			}
		}

		if (w1->start() >= w2->stop())
		{
			double Intersect_X, Intersect_Y;

			// judge the previous order of the two windows
			CalculateIntersectionPoint(w1->pseudo_x(), w1->pseudo_y(), w1->start(), 0, w2->pseudo_x(), w2->pseudo_y(), w2->stop(), 0, Intersect_X, Intersect_Y);

			face_pointer f = Tri.bottom_edge->opposite_face(Tri.face);
			edge_pointer e = f->next_edge(Tri.bottom_edge, Tri.left_vertex);
			double angle = f->vertex_angle(Tri.left_vertex);
			double Cx = e->length() * cos(angle);
			double Cy = e->length() * -sin(angle);

			if ((PointInTriangle(Intersect_X, Intersect_Y, Tri.bottom_edge->length(), Cx, Cy))
				&& (Intersect_Y <= 0) && (Intersect_Y >= w1->pseudo_y()) && (Intersect_Y >= w2->pseudo_y()))
			{
				double d1, d2;
				d1 = w1->d() + sqrt((Intersect_X - w1->pseudo_x()) * (Intersect_X - w1->pseudo_x()) + (Intersect_Y - w1->pseudo_y()) * (Intersect_Y - w1->pseudo_y()));
				d2 = w2->d() + sqrt((Intersect_X - w2->pseudo_x()) * (Intersect_X - w2->pseudo_x()) + (Intersect_Y - w2->pseudo_y()) * (Intersect_Y - w2->pseudo_y()));

				if (d1 < d2 * NUMERCIAL_EPSILON)
					return w2_invalid;
				if (d2 < d1 * NUMERCIAL_EPSILON)
					return w1_invalid;
			}
		}

		if (w2->start() >= w1->stop())
		{
			double Intersect_X, Intersect_Y;

			// judge the previous order of the two windows
			CalculateIntersectionPoint(w2->pseudo_x(), w2->pseudo_y(), w2->start(), 0, w1->pseudo_x(), w1->pseudo_y(), w1->stop(), 0, Intersect_X, Intersect_Y);

			face_pointer f = Tri.bottom_edge->opposite_face(Tri.face);
			edge_pointer e = f->next_edge(Tri.bottom_edge, Tri.left_vertex);
			double angle = f->vertex_angle(Tri.left_vertex);
			double Cx = e->length() * cos(angle);
			double Cy = e->length() * -sin(angle);

			if ((PointInTriangle(Intersect_X, Intersect_Y, Tri.bottom_edge->length(), Cx, Cy))
				&& (Intersect_Y <= 0) && (Intersect_Y >= w1->pseudo_y()) && (Intersect_Y >= w2->pseudo_y()))
			{
				double d1, d2;
				d1 = w1->d() + sqrt((Intersect_X - w1->pseudo_x()) * (Intersect_X - w1->pseudo_x()) + (Intersect_Y - w1->pseudo_y()) * (Intersect_Y - w1->pseudo_y()));
				d2 = w2->d() + sqrt((Intersect_X - w2->pseudo_x()) * (Intersect_X - w2->pseudo_x()) + (Intersect_Y - w2->pseudo_y()) * (Intersect_Y - w2->pseudo_y()));

				if (d1 < d2 - NUMERCIAL_EPSILON)
					return w2_invalid;
				if (d2 < d1 - NUMERCIAL_EPSILON)
					return w1_invalid;
			}
		}

		return both_valid;
	}

	inline void GeodesicAlgorithmExact::pairwise_windows_checking(list_pointer &list)
	{
		if (list->empty()) return;

		interval_pointer iter = list->begin();
		interval_pointer next, iter_t;

		next = iter->next();

		// traverse successive pairs of windows
		while ((!list->empty()) && (next != NULL))
		{
			windows_state ws = check_between_two_windows(iter, next);

			switch (ws)
			{
			case geodesic::w1_invalid:
				iter_t = iter;
				if (iter == list->begin())
				{
					iter = iter->next();
				}
				else
				{
					iter = iter->previous();
				}

				list->erase(iter_t);
				delete iter_t;
				--m_windows_wavefront;
				break;

			case geodesic::w2_invalid:
				list->erase(next);
				delete next;
				--m_windows_wavefront;
				break;

			case geodesic::both_valid:
				iter = iter->next();
				break;

			default:
				break;
			}

			next = iter->next();
		}
	}

	//------------------------- main operation ----------------------------
	inline void GeodesicAlgorithmExact::propagate_one_windows_list(list_pointer &list)
	{
		if (list->empty()) return;

		if (list->edge()->adjacent_faces().size() > 1)
		{
			// Rule 2: pairwise windows checking
			check_with_vertices(list);
			pairwise_windows_checking(list);

			// Rule 1: "One Angle Two Sides"
			find_separating_point(list);
			propagate_windows_to_two_edges(list);

		}
	}

	//-------------------------- main entry --------------------------
	inline void GeodesicAlgorithmExact::propagate(unsigned source)
	{
		// initialization
		m_source = source;
		initialize_propagation_data();

		clock_t start = clock();

		while (!m_vertex_queue.empty())
		{
			// (1) pop a vertex from M_VERTEX_QUEUE
			vertex_pointer vert = *m_vertex_queue.begin();
			m_vertex_queue.erase(m_vertex_queue.begin());

			// (2) update wavefront
			vert->state() = Vertex::INSIDE;
			for (unsigned i = 0; i < vert->adjacent_edges().size(); ++i)
			{
				vertex_pointer vert_it = vert->adjacent_edges()[i]->opposite_vertex(vert);
				if (vert_it->state() == Vertex::OUTSIDE) vert_it->state() = Vertex::FRONT;
			}

			// (3) handle saddle vertex
			if (vert->saddle_or_boundary()) create_pseudo_source_windows(vert, false);

			// (4) push window lists on the wavefront incident to v into M_LIST_QUEUE
			for (unsigned i = 0; i < vert->adjacent_edges().size(); ++i)
			{
				edge_pointer edge_it = vert->adjacent_edges()[i];
				if (!interval_list_0(edge_it)->empty()) m_list_queue.push(interval_list_0(edge_it));
				if (!interval_list_1(edge_it)->empty()) m_list_queue.push(interval_list_1(edge_it));
			}

			for (unsigned i = 0; i < vert->adjacent_faces().size(); ++i)
			{
				edge_pointer   edge_it = vert->adjacent_faces()[i]->opposite_edge(vert);
				vertex_pointer vert_it = (edge_it->adjacent_faces().size() < 2) ? NULL : edge_it->opposite_face(vert->adjacent_faces()[i])->opposite_vertex(edge_it);
				if (edge_it->adjacent_faces().size() < 2 || vert_it->state() != Vertex::OUTSIDE)
				{
					if (!interval_list_0(edge_it)->empty()) m_list_queue.push(interval_list_0(edge_it));
					if (!interval_list_1(edge_it)->empty()) m_list_queue.push(interval_list_1(edge_it));
				}
			}

			// (5) propagate window lists in a FIFO order
			while (!m_list_queue.empty())
			{
				// pop an list from M_LIST_QUEUE
				list_pointer list = m_list_queue.front();
				m_list_queue.pop();

				bool is_boundary = calculate_triangle_parameters(list, Tri);

				if (!is_boundary)
				{
					// propagate the window list using Rule 1 and 2
					wl_left.clear(); wl_right.clear();
					propagate_one_windows_list(list);

					// merge windows lists
					if (!wl_left.empty())
					{
						// in VTP, both "PrimeMerge" and "SecondMerge" connect window lists in an order-free way
						if (!Tri.left_list->empty())
						{
							Tri.left_list->begin()->previous() = wl_left.end();
							wl_left.end()->next() = Tri.left_list->begin();
							Tri.left_list->begin() = wl_left.begin();
						}
						else
						{
							Tri.left_list->begin() = wl_left.begin();
							Tri.left_list->end() = wl_left.end();
						}

						// push updated list into M_LIST_QUEUE
						if (((Tri.left_edge->v0()->state() == Vertex::INSIDE) || (Tri.left_edge->v1()->state() == Vertex::INSIDE)) && (!Tri.left_list->empty()))
							m_list_queue.push(Tri.left_list);
					}

					if (!wl_right.empty())
					{
						// in VTP, both "PrimeMerge" and "SecondMerge" connect window lists in an order-free way
						if (!Tri.right_list->empty())
						{
							Tri.right_list->end()->next() = wl_right.begin();
							wl_right.begin()->previous() = Tri.right_list->end();
							Tri.right_list->end() = wl_right.end();
						}
						else
						{
							Tri.right_list->begin() = wl_right.begin();
							Tri.right_list->end() = wl_right.end();
						}

						// push updated list into M_LIST_QUEUE
						if (((Tri.right_edge->v0()->state() == Vertex::INSIDE) || (Tri.right_edge->v1()->state() == Vertex::INSIDE)) && (!Tri.right_list->empty()))
							m_list_queue.push(Tri.right_list);
					}
				}

				list->clear();
			}

			// statistics
			if (m_vertex_queue.size() > m_queue_max_size)
				m_queue_max_size = m_vertex_queue.size();
		}

		clock_t stop = clock();
		m_time_consumed = (static_cast<double>(stop) - static_cast<double>(start)) / CLOCKS_PER_SEC;

	}

	//---------------------- print statistics --------------------------
	inline void GeodesicAlgorithmExact::print_statistics()
	{
		GeodesicAlgorithmBase::print_statistics();

		double memory = sizeof(Interval);

		//std::cout << std::endl;
		std::cout << "Peak number of intervals on wave-front " << m_windows_peak << std::endl;
		std::cout << "uses about " << memory * m_windows_peak / 1e6 << "MB of memory" << std::endl;
		std::cout << "total interval propagation number " << m_windows_propagation << std::endl;
		std::cout << "maximum interval queue size is " << m_queue_max_size << std::endl;
	}

}		//geodesic

#endif
