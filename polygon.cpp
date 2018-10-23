#include "polygon.h"
#include "utilities.h"
#include <limits>
#include <set>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <algorithm>

void Contour::boundingbox (Point& min, Point& max)
{
	min.x = min.y = numeric_limits<double>::max ();
	max.x = max.y = -numeric_limits<double>::max ();
	Contour::iterator i = begin();
	while (i != end()) {
		if (i->x < min.x)
			min.x = i->x;
		if (i->x > max.x)
			max.x = i->x;
		if (i->y < min.y)
			min.y = i->y;
		if (i->y > max.y)
			max.y = i->y;
		++i;
	}
}

bool Contour::counterclockwise ()
{
	if (_precomputedCC)
		return _CC;
	_precomputedCC = true;
	double area = 0.0;
	for (unsigned int c = 0; c < nvertices () - 1; c++)
		area += vertex (c).x * vertex (c+1).y - vertex (c+1).x *  vertex (c).y;
	area += vertex (nvertices ()-1).x * vertex (0).y - vertex (0).x *  vertex (nvertices ()-1).y;
	return _CC = area >= 0.0;
}

void Contour::move (double x, double y)
{
	for (unsigned int i = 0; i < points.size (); i++) {
		points[i].x += x;
		points[i].y += y;
	}
}

ostream& operator<< (ostream& o, Contour& c)
{
	o << c.nvertices () << endl;
	Contour::iterator i = c.begin();
	while (i != c.end()) {
		o << '\t' << i->x << " " << i->y << endl;
		++i;
	}
	return o;
}

Polygon::Polygon (const string& filename)
{
	ifstream f (filename.c_str ());
	if (!f) {
		cerr << "Error opening " << filename << '\n';
		exit (1);
	}
	f >> *this;
	if (! f.eof ())
		std::cerr << "An error reading file " << filename << " happened\n";
}

unsigned Polygon::nvertices () const
{
	unsigned int nv = 0;
	for (unsigned int i = 0; i < ncontours (); i++)
		nv += contours[i].nvertices ();
	return nv;
}

void Polygon::boundingbox (Point& min, Point& max)
{
	min.x = min.y = numeric_limits<double>::max ();
	max.x = max.y = -numeric_limits<double>::max ();
	Point mintmp;
	Point maxtmp;
	for (unsigned int i = 0; i < ncontours (); i++) {
		contours[i].boundingbox (mintmp, maxtmp);
		if (mintmp.x < min.x)
			min.x = mintmp.x;
		if (maxtmp.x > max.x)
			max.x = maxtmp.x;
		if (mintmp.y < min.y)
			min.y = mintmp.y;
		if (maxtmp.y > max.y)
			max.y = maxtmp.y;
	}
}

void Polygon::move (double x, double y)
{
	for (unsigned int i = 0; i < contours.size (); i++)
		contours[i].move (x, y);
}

ostream& operator<< (ostream& o, Polygon& p)
{
	o << p.ncontours () << std::endl;
	for (unsigned int i = 0; i < p.ncontours (); i++)   // write the contours
		o << p.contour (i);
	for (unsigned int i = 0; i < p.ncontours (); i++) { // write the holes of every contour
		if (p.contour (i).nholes () > 0) {
			o << i << ": ";
			for (unsigned int j = 0; j < p.contour (i).nholes (); j++)
				o << p.contour (i).hole (j) << (j == p.contour (i).nholes () - 1 ? '\n' : ' ');
		}
	}
	return o;
}

istream& operator>> (istream& is, Polygon& p)
{
	// read the contours
	int ncontours;
	double px, py;
	is >> ncontours;
	for (int i = 0; i < ncontours; i++) {
		int npoints;
		is >> npoints;
		p.push_back (Contour ());
		Contour& contour = p.back ();
		for (int j = 0; j < npoints; j++) {
			is >> px >> py;
			if (j > 0 && px == contour.back ().x && py == contour.back ().y)
				continue;
			if (j == npoints-1 && px == contour.vertex (0).x && py == contour.vertex (0).y)
				continue;
			contour.add (Point (px, py));
		}
		if (contour.nvertices () < 3) {
			p.pop_back ();
			continue;
		}
	}
	// read holes information
	int contourId;
	char aux;
	std::string restOfLine;
	while (is >> contourId) {
		is >> aux; // read the character :
		if (aux != ':')
			break;
		std::getline (is, restOfLine);
		std::istringstream iss (restOfLine);
		int hole;
		while (iss >> hole) {
			p[contourId].addHole (hole);
			p[hole].setExternal (false);
		}
		if (! iss.eof ())
			break;
	}
	return is;
}

/*
 * The following code is necessary for implementing the computeHoles member function
 * 
 */

namespace { // start of anonymous namespace
	struct SweepEvent;
	struct SegmentComp : public binary_function<SweepEvent*, SweepEvent*, bool> {
		bool operator() (SweepEvent* e1, SweepEvent* e2);
	};

	struct SweepEvent {
		Point p;   // point associated with the event
		bool left; // is the point the left endpoint of the segment (p, other->p)?
		int pl;    // Polygon to which the associated segment belongs to
		SweepEvent* other; // Event associated to the other endpoint of the segment
		/**  Does the segment (p, other->p) represent an inside-outside transition in the polygon for a vertical ray from (p.x, -infinite) that crosses the segment? */
		bool inOut;
		set<SweepEvent*, SegmentComp>::iterator poss; // Only used in "left" events. Position of the event (segment) in S

		/** Class constructor */
		SweepEvent (const Point& pp, bool b, int apl) : p (pp), left (b), pl (apl) {}
		/** Return the segment associated to the SweepEvent */
		Segment segment () { return Segment (p, other->p); }
		/** Is the segment (p, other->p) below point x */
		bool below (const Point& x) const { return (left) ? signedArea (p, other->p, x) > 0 : signedArea (other->p, p, x) > 0; }
		/** Is the segment (p, other->p) above point x */
		bool above (const Point& x) const { return !below (x); }
	};

	struct SweepEventComp : public binary_function<SweepEvent*, SweepEvent*, bool> {
		bool operator() (SweepEvent* e1, SweepEvent* e2) {
			if (e1->p.x < e2->p.x) // Different x coordinate
				return true;
			if (e2->p.x < e1->p.x) // Different x coordinate
				return false;
			if (e1->p != e2->p) // Different points, but same x coordinate. The event with lower y coordinate is processed first
				return e1->p.y < e2->p.y;
			if (e1->left != e2->left) // Same point, but one is a left endpoint and the other a right endpoint. The right endpoint is processed first
				return !e1->left;
			// Same point, both events are left endpoints or both are right endpoints. The event associate to the bottom segment is processed first
			return e1->below (e2->other->p);
		}
};

} // end of anonymous namespace

ostream& operator<< (ostream& o, SweepEvent& e)
{
	return o << " Point: " << e.p << " Other point: " << e.other->p << (e.left ? " (Left) " : " (Right) ")
            << (e.inOut ? " (In-Out) " : " (Out-In) ") << " Polygon: " << e.pl;
}

bool SegmentComp::operator() (SweepEvent* e1, SweepEvent* e2) {
	if (e1 == e2)
		return false;
	if (signedArea (e1->p, e1->other->p, e2->p) != 0 || signedArea (e1->p, e1->other->p, e2->other->p) != 0) {
		// Segments are not collinear
		// If they share their left endpoint use the right endpoint to sort
		if (e1->p == e2->p)
			return e1->below (e2->other->p);
		// Different points
		SweepEventComp comp;
		if (comp (e1, e2))  // has the segment associated to e1 been sorted in evp before the segment associated to e2?
			return e1->below (e2->p);
		// The segment associated to e2 has been sorted in evp before the segment associated to e1
		return e2->above (e1->p);
	} 
	// Segments are collinear. Just a consistent criterion is used
	if (e1->p == e2->p)
		return e1 < e2;
	SweepEventComp comp;
	return comp (e1, e2);
}

//#define _DEBUG_
void Polygon::computeHoles ()
{
	if (ncontours () < 2) {
		if (ncontours () == 1 && contour (0).clockwise ())
			contour (0).changeOrientation ();
		return;
	}
	vector<SweepEvent> ev;
	vector<SweepEvent*> evp;
	ev.reserve (nvertices ()*2);
	evp.reserve (nvertices ()*2);
	for (unsigned i = 0; i < ncontours (); i++) {
//		cout << contour (i);
		contour (i).setCounterClockwise ();
//		cout << contour (i);
		for (unsigned j = 0; j < contour (i).nedges (); j++) {
			Segment s = contour(i).segment (j);
			if (s.begin ().x == s.end ().x) // vertical segments are not processed
				continue;
			ev.push_back (SweepEvent (s.begin (), true, i));
			ev.push_back (SweepEvent (s.end (), true, i));
			SweepEvent* se1 = &ev[ev.size ()-2];
			SweepEvent* se2 = &ev[ev.size ()-1];
			se1->other = se2;
			se2->other = se1;
			if (se1->p.x < se2->p.x) {
				se2->left = false;
				se1->inOut = false;
			} else {
				se1->left = false;
				se2->inOut = true;
			}
			evp.push_back (se1);
			evp.push_back (se2);
		}
	}
	sort (evp.begin (), evp.end (), SweepEventComp ()); 

	set<SweepEvent*, SegmentComp> S; // Status line
	vector<bool> processed (ncontours (), false);
	vector<int> holeOf (ncontours (), -1);
	int nprocessed = 0;
	for (unsigned i = 0; i < evp.size () && nprocessed < ncontours (); i++)  {
		SweepEvent* e = evp[i];
		#ifdef _DEBUG_
		cout << "Process event: " << *e << endl;
		#endif

		if (e->left) { // the segment must be inserted into S
			e->poss = S.insert(e).first;
			if (!processed[e->pl]) {
				processed[e->pl] = true;
				nprocessed++;
				set<SweepEvent*, SegmentComp>::iterator prev = e->poss;
				if (prev == S.begin()) {
					contour (e->pl).setCounterClockwise ();
				} else {
					prev--;
					if (!(*prev)->inOut) {
						holeOf[e->pl] = (*prev)->pl;
						contour (e->pl).setExternal (false);
						contour ((*prev)->pl).addHole (e->pl);
						if (contour((*prev)->pl).counterclockwise ())
							contour (e->pl).setClockwise ();
						else
							contour (e->pl).setCounterClockwise ();
					} else if (holeOf[(*prev)->pl] != -1) {
						holeOf[e->pl] = holeOf[(*prev)->pl];
						contour (e->pl).setExternal (false);
						contour (holeOf[e->pl]).addHole (e->pl);
						if (contour(holeOf[e->pl]).counterclockwise ())
							contour (e->pl).setClockwise ();
						else
							contour (e->pl).setCounterClockwise ();
					} else {
						contour (e->pl).setCounterClockwise ();
					}
				}
			}
		} else { // the segment must be removed from S
			S.erase (e->other->poss);
		}
		#ifdef _DEBUG_
		cout << "Tras ajuste: " << endl;
		for (set<SweepEvent*, SegmentComp>::const_iterator it2 = S.begin(); it2 != S.end(); it2++)
			cout << **it2 << endl;
		cout << endl;
		string st;
		getline (cin, st);
		#endif
	}
}
