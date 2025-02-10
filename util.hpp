#ifndef UTILS
#define UTILS


#include <cmath>
#include <vector>
#include <cstdio>
#include <random>
#include <string>
#include <cstring>
#include <fstream>
#include <iostream>
#include <algorithm>

namespace utils{
class Point {
public:

    Point (){}
    Point(double x, double y);

    Point operator+ (Point& first);

    Point operator- (Point& first);

    Point operator+= (Point& first);

    Point operator-= (Point& first);

    Point operator-();

    double length();

    double getX() const;

    double getY() const;

    bool operator< (const Point& other) const {
        return other.x < this->x || (other.x == this->x && other.y < this->y);
    }

    double x,y;

};

class Line {

public:
    Line(double m, double b);

    static Line buildByPoints(Point& start, Point& end);

    static Line buildByPointAndAngle(Point& start, double angle);

    double getM() const;

    double getB() const;

private:
    double m, b;
};


class PointUtil {

public:
    static double orientation(Point& one, Point& two, Point& three);

    static Point vector(double angle, double length);

    static Point perpendicular(Point& one, Point& two, double length, int orientation);

    static const int CLOCKWISE = 1;
    static const int COUNTERCLOCKWISE = -1;
};


class LineSegment {

public:
    LineSegment(Point &start, Point &end);

    LineSegment(const Line &line, const Point &start, const Point &end);

    LineSegment(const LineSegment& copySegment);

    ~LineSegment();

    double length();

    Line getLine();

    Point getStart();

    Point* getStartPtr();

    Point getEnd();

    Point* getEndPtr();

    Line line;
    Point *start, *end;
};


class Ellipse {

public:
    Ellipse (){}
    Ellipse(Point centergiven, double x_radius, double y_radius)
    {
        center = centergiven;
        radius_x = x_radius;
        radius_y = y_radius;
        _size = M_PI * radius_x * radius_y;
    }

    bool inside(const Point &vector);

    double size();

    LineSegment segmentIntersections(LineSegment &segment);

    LineSegment intersections(Line &line);

    bool crosses(LineSegment &segment);

    bool crossesEdge(LineSegment &segment);

    Point getCross(LineSegment &segment) ;

    double edgeGradient(Point& point);

    Point getCenter() {
        return center;
    }

    double getXRadius() {
        return radius_x;
    }

    double getYRadius() {
        return radius_y;
    }
    
    Point center;
    double radius_x, radius_y;
    double _size;
};

double get_dist (Point A, Point B)
{
    return sqrt ((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y));
}

double gaussian (Point &input, Point &gaussianCenter, Point &gaussianVar)
{
    double X = input.x - gaussianCenter.x;
    double Y = input.y - gaussianCenter.y;
    double varX = gaussianVar.x;
    double varY = gaussianVar.y;
    return exp ( - (X * X / (2*varX)) - (Y * Y/(2*varY)));
}

std::vector<double> getGaussian (std::vector<Point> &points,
                            std::vector<Point> &gaussianCenters,
                            std::vector<Point> &gaussianVars)
{
    std::vector<double> fval;
    
    for (int i = 0; i < points.size(); i ++)
    {
        double sum = 0;
        for (int j = 0; j < gaussianCenters.size(); j++)
            sum = sum + gaussian(points[i], gaussianCenters[j], gaussianVars[j]);
        fval.push_back (sum);
    }
    
    return fval;
}


std::vector<double> get_Gaussian_vector (std::vector<Point> points,
                                    std::vector<Point> &gaussianCenters, 
                                    std::vector<Point> &gaussianVars,
                                    int id)
{
    std::vector<double> x;
    std::vector<double> fval = getGaussian(points, gaussianCenters, gaussianVars);
    //below is buggy, does not incorporate variance
    x.push_back (-2*(points[1].x-gaussianCenters[id].x) * fval[1]);
    x.push_back (-2*(points[1].y-gaussianCenters[id].y) * fval[1]);
    return x;
}

Point::Point(double x, double y) : x(x), y(y) {}

Point Point::operator+(Point &first) {
    return Point(x + first.x, y + first.y);
}

Point Point::operator+=(Point &first){
    this->x += first.x;
    this->y += first.y;
    return *this;
}

Point Point::operator-(Point &first) {
    return Point(x - first.x, y - first.y);
}

Point Point::operator-=(Point &first){
    this->x -= first.x;
    this->y -= first.y;
    return *this;
}

Point Point::operator-() {
    return Point(-this->x, -this->y);
}

double Point::length() {
    return sqrt(x * x + y * y);
}

double Point::getX() const {
    return x;
}

double Point::getY() const {
    return y;
}

double PointUtil::orientation(Point& one, Point& two, Point& three) {
    double k=(two.getY() - one.getY())*(three.getX() - two.getX())-(two.getX() - one.getX()) * (three.getY() - two.getY());

    if(k>0) {
        return CLOCKWISE;
    } else {

        return COUNTERCLOCKWISE;
    }
}

Point PointUtil::vector(double angle, double length) {
    return Point(length * cos(angle), length * sin(angle));
}

Point PointUtil::perpendicular(Point &one, Point &two, double length, int orientation) {
    double delta_x = two.getX() - one.getX();
    double delta_y = two.getY() - one.getY();
    double angle = atan2(delta_y, delta_x);
    return vector(angle + (orientation * M_PI_2), length);
}


Line::Line(double m, double b) : m(m), b(b) {}

Line Line::buildByPoints(Point &start, Point &end) {
    double m = (end.getY() - start.getY()) / (end.getX() - start.getX() + 1e-9); //divide by zero case solved by 1e-9
    double b = start.getY() - (m * start.getX());

    return Line(m, b);
}

Line Line::buildByPointAndAngle(Point &start, double angle) {
    double m = tan(angle);
    double b = start.getY() - (m * start.getX());

    return Line(m, b);
}

double Line::getM() const {
    return m;
}

double Line::getB() const {
    return b;
}

LineSegment::LineSegment(Point &start, Point &end) : line(Line::buildByPoints(start, end)), start(new Point(start.getX(), start.getY())), end(new Point(end.getX(), end.getY())) {}

LineSegment::LineSegment(const Line &line, const Point &start, const Point &end) : line(line), start(new Point(start.getX(), start.getY())), end(new Point(end.getX(), end.getY())) {}

LineSegment::LineSegment(const LineSegment &copySegment): line(copySegment.line), start(new Point(copySegment.start->getX(), copySegment.start->getY())), end(new Point(copySegment.end->getX(), copySegment.end->getY())) {}

double LineSegment::length() {
    Point vector = (*end - *start);
    return vector.length();
}

Line LineSegment::getLine() {
    return line;
}

Point LineSegment::getStart() {
    return *start;
}

Point* LineSegment::getStartPtr() {
    return start;
}

Point LineSegment::getEnd() {
    return *end;
}

Point* LineSegment::getEndPtr() {
    return end;
}

LineSegment::~LineSegment() {
    delete start;
    delete end;
}

bool Ellipse::inside(const Point &vector) {
    return pow(((vector.getX() - center.getX()) / radius_x), 2) +
           pow(((vector.getY() - center.getY()) / radius_y), 2) <= 1;
}

double Ellipse::size() {
    return _size;
}

LineSegment Ellipse::segmentIntersections(LineSegment &segment) {
    Line line = segment.getLine();
    LineSegment intersectionSegment = intersections(line);

    Point start = intersectionSegment.getStart();
    Point end = intersectionSegment.getEnd();

    if(start.getX() > end.getX()) {
        Point tmp = start;
        start = end;
        end = tmp;
    }

    Point segmentStart = segment.getStart();
    Point segmentEnd = segment.getEnd();

    if(segmentStart.getX() > segmentEnd.getX()) {
        Point tmp = segmentStart;
        segmentStart = segmentEnd;
        segmentEnd = tmp;
    }

    if(segmentStart.getX() > start.getX()) {
        start = segmentStart;
    }
    if(segmentEnd.getX() < end.getX()) {
        end = segmentEnd;
    }
    if(start.getX() > end.getX()) {
        start = end;
    }

    return LineSegment(segment.getLine(), start, end);
}

LineSegment Ellipse::intersections(Line &line) {
    double rx = radius_x * radius_x;
    double ry = radius_y * radius_y;

    double a = (1 / rx) + (line.getM() * line.getM() / ry);
    double b = (2 * line.getB() * line.getM() / ry) - (2 * center.getX() / rx) - (2 * center.getY() * line.getM() / ry);
    double c = (line.getB() * line.getB() / ry) - (2 * line.getB() * center.getY() / ry) + (center.getX() * center.getX() / rx) + (center.getY() * center.getY() / ry) - 1;

    // Solution using Quadratic equation -b +- sqrt(b^2 - 4ac)/2a
    // where ax^2 + bx + c = 0
    double discriminant = pow(b,2) - (4 * a * c);

    if (discriminant > 0){
        double x1 = ((-b) + sqrt(discriminant)) / (2 * a);
        double x2 = ((-b) - sqrt(discriminant)) / (2 * a);

        double y1 = (line.getM() * x1) + line.getB();
        double y2 = (line.getM() * x2) + line.getB();

        return LineSegment(line, Point(x1, y1), Point(x2, y2));
    } else {
        return LineSegment(line, Point(0, 0), Point(0, 0));
    }
}

bool Ellipse::crosses(LineSegment &segment) {
    Line line = segment.getLine();

    double rx = radius_x * radius_x;
    double ry = radius_y * radius_y;

    double a = (1 / rx) + (line.getM() * line.getM() / ry);
    double b = (2 * line.getB() * line.getM() / ry) - (2 * center.getX() / rx) - (2 * center.getY() * line.getM() / ry);
    double c = (line.getB() * line.getB() / ry) - (2 * line.getB() * center.getY() / ry) + (center.getX() * center.getX() / rx) + (center.getY() * center.getY() / ry) - 1;

    // Solution using Quadratic equation -b +- sqrt(b^2 - 4ac)/2a
    // where ax^2 + bx + c = 0
    double discriminant = pow(b,2) - (4 * a * c);

    if (discriminant > 0){
        double x1 = ((-b) + sqrt(discriminant)) / (2 * a);
        double x2 = ((-b) - sqrt(discriminant)) / (2 * a);

        return (segment.getStart().getX() < x1 && segment.getStart().getX() > x2)
            || (segment.getEnd().getX() < x1 && segment.getEnd().getX() > x2);
    } else {
        return false;
    }
}


bool Ellipse::crossesEdge(LineSegment &segment) {
    Line line = segment.getLine();
    LineSegment intersectionSegment = intersections(line);

    Point start = intersectionSegment.getStart();
    Point end = intersectionSegment.getEnd();

    Point segmentStart = segment.getStart();
    Point segmentEnd = segment.getEnd();
    
    //check if discriminant is negative
    if (start.x == 0 && start.y == 0 && end.x == 0 && end.y == 0)
        return false;

    if(segmentStart.getX() > segmentEnd.getX()) {
        Point tmp = segmentStart;
        segmentStart = segmentEnd;
        segmentEnd = tmp;
    }

    return (start.getX() > segmentStart.getX() && start.getX() < segmentEnd.getX()) ||
            (end.getX() > segmentStart.getX() && end.getX() < segmentEnd.getX());
}


double Ellipse::edgeGradient(Point& point) {
    return atan2(-pow(radius_y, 2) * (point.getX()-center.getX()), (pow(radius_x, 2) * (point.getY()-center.getY())));
}
}

#endif