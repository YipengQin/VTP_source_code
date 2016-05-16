#ifndef GEODESIC_MESH_ELEMENTS
#define GEODESIC_MESH_ELEMENTS

// here we define the building elements of the mesh: 
// 3D-points, vertices, edges, faces, and surface points

#include "stdafx.h"
#include "geodesic_constants_and_simple_functions.h"

namespace geodesic{

class Vertex;
class Edge;
class Face;
class Mesh;
class MeshElementBase;

typedef Vertex* vertex_pointer;
typedef Edge* edge_pointer;
typedef Face* face_pointer;
typedef Mesh* mesh_pointer;
typedef MeshElementBase* base_pointer;

// Base Elements

template <class Data>		//simple vector that stores info about mesh references 
class SimpleVector			//for efficiency, it uses an outside memory allocator
{
public:
	SimpleVector():
	  m_size(0),
	  m_begin(NULL)
	{};

	typedef Data* iterator;

	unsigned size(){return m_size;};
	iterator begin(){return m_begin;};
	iterator end(){return m_begin + m_size;};

	template<class DataPointer>
	void set_allocation(DataPointer begin, unsigned size)
	{
		assert(begin != NULL || size == 0);
		m_size = size;
		m_begin = (iterator)begin;
	}

	Data& operator[](unsigned i)
	{
		assert(i < m_size);
		return *(m_begin + i);
	}

	void clear()
	{
		m_size = 0;
		m_begin = NULL;
	}

private:
	unsigned m_size;
	Data* m_begin;
};

enum PointType
{
    VERTEX,
    EDGE,
    FACE,
	UNDEFINED_POINT
};

class MeshElementBase	//prototype of vertices, edges and faces
{
public:
	typedef SimpleVector<vertex_pointer> vertex_pointer_vector;	
	typedef SimpleVector<edge_pointer>   edge_pointer_vector;	
	typedef SimpleVector<face_pointer>   face_pointer_vector;	

	MeshElementBase():
		m_id(0),
		m_type(UNDEFINED_POINT)
	{};

	vertex_pointer_vector& adjacent_vertices(){return m_adjacent_vertices;};
	edge_pointer_vector& adjacent_edges(){return m_adjacent_edges;};
	face_pointer_vector& adjacent_faces(){return m_adjacent_faces;};

	unsigned& id(){return m_id;}; 
	PointType type(){return m_type;};

protected:
	vertex_pointer_vector m_adjacent_vertices;		//list of the adjacent vertices
	edge_pointer_vector m_adjacent_edges;			//list of the adjacent edges
	face_pointer_vector m_adjacent_faces;			//list of the adjacent faces

	unsigned m_id;							//unique id
	PointType m_type;							//vertex, edge or face
};

class Point3D			//point in 3D and corresponding operations
{
public:
	Point3D(){};
	Point3D(Point3D* p)
	{
		x() = p->x();
		y() = p->y();
		z() = p->z();
	};

	double* xyz(){return m_coordinates;};
	double& x(){return *m_coordinates;};
	double& y(){return *(m_coordinates+1);};
	double& z(){return *(m_coordinates+2);};

	void set(double new_x, double new_y, double new_z)
	{
		x() = new_x;
		y() = new_y;
		z() = new_z;
	}

	void set(double* data)
	{
		x() = *data;
		y() = *(data+1);
		z() = *(data+2);
	}

	double distance(double* v)
	{
		double dx = m_coordinates[0] - v[0];
		double dy = m_coordinates[1] - v[1];
		double dz = m_coordinates[2] - v[2];

		return sqrt(dx*dx + dy*dy + dz*dz);
	};

    double distance(Point3D* v)
	{
		return distance(v->xyz());
	};

	void add(Point3D* v)
	{
		x() += v->x();
		y() += v->y();
		z() += v->z();
	};

	void multiply(double v)
	{
		x() *= v;
		y() *= v;
		z() *= v;
	};

private:
	double m_coordinates[3];					//xyz
};

// Mesh Elements

class Vertex: public MeshElementBase, public Point3D
{

public:

	enum VertexState
	{
		OUTSIDE,
		FRONT,
		INSIDE
	};

	Vertex()
	{
		m_type = VERTEX;
		m_geodesic_distance = geodesic::GEODESIC_INF;  // Initialize Distance
		m_state = OUTSIDE;
		m_incident_face = NULL;
		m_saddle_or_boundary = false;
		m_incident_point = NULL; 
	};

	~Vertex(){};

	bool operator()(vertex_pointer const &x, vertex_pointer const &y) const // compare geodesic distances of two vertices
	{
		return (x->geodesic_distance() < y->geodesic_distance());
	}

	// member operations
	bool&      saddle_or_boundary() { return m_saddle_or_boundary; }; // return whether the vertex is saddle vertex or boundary vertex

	face_pointer&   incident_face() { return m_incident_face; }; // return the incident face of shortest path to the vertex
	double&        incident_point() { return m_incident_point; }; // return the incident point position of the shortest path
	
	VertexState&            state() { return m_state; }; // return the state of vertex according to traversed area R
	double&     geodesic_distance() { return m_geodesic_distance; }; // return the geodesic distance of the vertex

private:	

	// member
	bool           m_saddle_or_boundary;		//it is true if total adjacent angle is larger than 2*PI or this vertex belongs to the mesh boundary

	face_pointer   m_incident_face; // the incident face of shortest path to the vertex
	double         m_incident_point; // the incident point position of the shortest path

	VertexState    m_state; // the state of the vertex
	double         m_geodesic_distance; // the geodesic distance of the vertex
		
};

class Face: public MeshElementBase
{
public:
	Face()
	{
		m_type = FACE;
	};

	~Face(){};

	vertex_pointer opposite_vertex(edge_pointer e);
	edge_pointer   opposite_edge(vertex_pointer v);
	edge_pointer   next_edge(edge_pointer e, vertex_pointer v);

	double vertex_angle(vertex_pointer v)
	{
		for(unsigned i = 0; i < 3; ++i)
		{
			if(adjacent_vertices()[i]->id() == v->id())
			{
				return m_corner_angles[i];
			}
		}
		assert(0);
		return 0;
	}

	double* corner_angles(){return m_corner_angles;};

private:
	double m_corner_angles[3];		//triangle angles in radians; angles correspond to vertices in m_adjacent_vertices
};

class Edge: public MeshElementBase
{
public:
	Edge()
	{
		m_type = EDGE;
	};

	~Edge(){};

	double& length(){return m_length;};

	face_pointer opposite_face(face_pointer f)
	{
		if(adjacent_faces().size() == 1)
		{
			assert(adjacent_faces()[0]->id() == f->id());
			return NULL;
		}

		assert(adjacent_faces()[0]->id() == f->id() || 
			   adjacent_faces()[1]->id() == f->id());

		return adjacent_faces()[0]->id() == f->id() ? 
			   adjacent_faces()[1] : adjacent_faces()[0];
	};

	vertex_pointer opposite_vertex(vertex_pointer v)
	{
		assert(belongs(v));

		return adjacent_vertices()[0]->id() == v->id() ?
			   adjacent_vertices()[1] : adjacent_vertices()[0];
	};

	bool belongs(vertex_pointer v)
	{
		return adjacent_vertices()[0]->id() == v->id() || 
			   adjacent_vertices()[1]->id() == v->id();
	}

	bool is_boundary(){return adjacent_faces().size() == 1;};

	vertex_pointer v0(){return adjacent_vertices()[0];};
	vertex_pointer v1(){return adjacent_vertices()[1];};

	void local_coordinates(Point3D* point, 
						   double& x, 
						   double& y)
	{
		double d0 = point->distance(v0());
		if(d0 < 1e-50)
		{
			x = 0.0;
			y = 0.0;
			return;
		}

		double d1 = point->distance(v1());
		if(d1 < 1e-50)
		{
			x = m_length;
			y = 0.0;
			return;
		}

		x = m_length/2.0 + (d0*d0 - d1*d1)/(2.0*m_length);
		y = sqrt(std::max(0.0, d0*d0 - x*x));
		return;
	}

private:
	double m_length;							//length of the edge
};

class SurfacePoint:public Point3D  //point on the surface of the mesh
{
public:
	SurfacePoint():
		m_p(NULL)
	{};

	SurfacePoint(vertex_pointer v):		//set the surface point in the vertex
		SurfacePoint::Point3D(v),
		m_p(v)
	{};

	SurfacePoint(face_pointer f):		//set the surface point in the center of the face
		m_p(f)
	{
		set(0,0,0);
		add(f->adjacent_vertices()[0]);
		add(f->adjacent_vertices()[1]);
		add(f->adjacent_vertices()[2]);
		multiply(1./3.);
	};

	SurfacePoint(edge_pointer e,		//set the surface point in the middle of the edge
				 double a = 0.5):		
		m_p(e)
	{
		double b = 1 - a;

		vertex_pointer v0 = e->adjacent_vertices()[0];
		vertex_pointer v1 = e->adjacent_vertices()[1];

		x() = b*v0->x() + a*v1->x();
		y() = b*v0->y() + a*v1->y();
		z() = b*v0->z() + a*v1->z();
	};

	SurfacePoint(base_pointer g, 
				 double x,
				 double y,
				 double z,
				 PointType t = UNDEFINED_POINT):
		m_p(g)
	{
		set(x,y,z);
	};

	void initialize(SurfacePoint const& p)
	{
		*this = p;
	}

	~SurfacePoint(){};

	PointType type(){return m_p ? m_p->type() : UNDEFINED_POINT;};
	base_pointer& base_element(){return m_p;};
protected:
	base_pointer m_p;			//could be face, vertex or edge pointer
};

// Face Member Functions

inline edge_pointer Face::opposite_edge(vertex_pointer v)
{
	for(unsigned i = 0; i < 3; ++i)
	{
		edge_pointer e = adjacent_edges()[i];
		if(!e->belongs(v))
		{
			return e;
		}
	}
	assert(0);
	return NULL;
}

inline vertex_pointer Face::opposite_vertex(edge_pointer e)
{
	for(unsigned i = 0; i < 3; ++i)
	{
		vertex_pointer v = adjacent_vertices()[i];
		if(!e->belongs(v))
		{
			return v;
		}
	}
	assert(0);
	return NULL;
}

inline edge_pointer Face::next_edge(edge_pointer e, vertex_pointer v)
{
	assert(e->belongs(v));

	for(unsigned i = 0; i < 3; ++i)
	{
		edge_pointer next = adjacent_edges()[i];
		if(e->id() != next->id() && next->belongs(v))
		{
			return next;
		}
	}
	assert(0);
	return NULL;
}

// HalfEdge operations

struct HalfEdge			//prototype of the edge; used for mesh construction
{
	unsigned face_id;
	unsigned vertex_0;		//adjacent vertices sorted by id value
	unsigned vertex_1;		//they are sorted, vertex_0 < vertex_1
};

inline bool operator < (const HalfEdge &x, const HalfEdge &y)
{
	if(x.vertex_0 == y.vertex_0)
	{
	    return x.vertex_1 < y.vertex_1;
	}
	else
	{
		return x.vertex_0 < y.vertex_0;
	}
}

inline bool operator != (const HalfEdge &x, const HalfEdge &y)
{
	return x.vertex_0 != y.vertex_0 || x.vertex_1 != y.vertex_1;
}

inline bool operator == (const HalfEdge &x, const HalfEdge &y)
{
	return x.vertex_0 == y.vertex_0 && x.vertex_1 == y.vertex_1;
}

} //geodesic

#endif	
