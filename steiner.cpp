#include <iostream>
#include <fstream>
#include <cstdint>
#include <list>
#include <array>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <algorithm>
#include <functional>
#include <array>
#include <iterator>
#include <string>
#include <regex>
#include <unordered_map>
#include <set>
#include <queue>
#include <utility>
#include <algorithm>


enum Type
{
    NONE,
    TERM,
    PATH
};

struct Point
{
    int x;
    int y;
    
    int status;
    bool operator==(const Point& that) 
    {
        if(this->x == that.x && this->y == that.y)
            return true;
        return false;
    }

    bool operator<(const Point& that) const
    {
        if(this->y < that.y)
            return true;
       else if(this->y == that.y)
       {
            if(this->x < that.x)
                return true;
       }    
       return false;
            
    }

    friend std::ostream& operator<<(std::ostream& os, const Point& point);
};

std::ostream& operator<<(std::ostream& os, const Point& point)
{
    os <<  "(" << point.x << ", " << point.y << ")";
    return os;
}

struct Edge
{
    int u;
    int v;
    int weight;

};

std::vector<Point> parsePoints(std::ifstream& file)
{
    std::vector<Point> points;
    std::string temp;
    std::regex re(".*point.*x=\"([0-9]+).*y=\"([0-9]+)\".*");

    while(std::getline(file, temp))
    {

        std::smatch match;
        if(std::regex_search(temp, match, re))
        {
            int x = std::stoi(match[1]);
            int y = std::stoi(match[2]);
            Point p = {x, y};
            points.push_back(p);
        }
   
    }
    
    return points;
}

class Graph
{
public:
    //Graph
    int addNode(Point point);
    void addEdge(int i, int j);
    void dump();
    int find(Point point);
    std::vector<int> getPointsId(std::vector<Point>& points);

    std::vector<Point> nodes;
    std::vector< std::vector<int> > adj;
    std::vector<int> terminals;
    std::set<Point> dupl;
    int size;
};

int Graph::find(Point point)
{
    for(int i = 0; i < nodes.size(); i++)
        if(nodes[i] == point)
            return i;
    return -1;
}

std::vector<int> Graph::getPointsId(std::vector<Point>& points)
{
    std::vector<int> ids;
    for(auto& p : points)
    {
        for(int i = 0; i < nodes.size(); i++)
        {
            if(p == nodes[i])
                ids.push_back(i);
        }
    }
    return ids;
}

int Graph::addNode(Point point)
{
    if(dupl.find(point) == dupl.end())
    {
        nodes.push_back(point);
        adj.push_back(std::vector<int>());
        size++;
        dupl.insert(point);
        return nodes.size() - 1;
    }
    else
    {
        for(int i = 0; i < nodes.size(); i++)
            if(nodes[i] == point)
                return i;
    }
}

void Graph::addEdge(int i, int j)
{
    adj[i].push_back(j);
    adj[j].push_back(i);
}

void Graph::dump()
{
    for(auto n : nodes)
        std::cout << n << ": " << n.status << "\n";
    std::cout << "\n";
    
    std::set< std::pair<int, int> > used;

    for(int i = 0; i < nodes.size(); i++)
    {
        for(auto& n : adj[i])
        {
            if(used.find({i, n}) == used.end())
                std::cout << nodes[i] << "--" << nodes[n] << "\n";
            used.insert({i, n});
            used.insert({n, i});
        }
    }

}

Graph* createGridGraph(int size)
{
    Graph* grid = new Graph();
    
    for(int i = 0; i < size; i++)
    {
        for(int j = 0; j < size; j++)
            grid->addNode({j, i, NONE}); 
    }
    
    for(int i = 0; i < size; i++)
    {
        for(int j = 0; j < size; j++)
        {
            if(i == 0)
            {
                if(j != size - 1)
                    grid->addEdge(i+j, i+j+1);
                grid->addEdge(i+j, i+j+size);

            }

            else if(i == size - 1)
            {
                if(j != size - 1)
                    grid->addEdge(i*size+j, i*size+j+1);
            }

            else
            {
                if(j != size - 1)
                    grid->addEdge(i*size+j, i*size+j+1);
                grid->addEdge(i*size+j, i*size+j+size);

            }
            
        }
    }
    
    return grid;
}

std::vector<int> addTerminal(Graph* grid, int src)
{
    std::queue<int> queue;
    assert(grid->size == grid->nodes.size());
    std::vector<bool> used(grid->size);
    std::vector<int> parents(grid->size);


    queue.push(src);
    used[src] = true;
    parents[src] = -1;

    int dst = -1;
    while(!queue.empty())
    {
        int cur = queue.front();
        queue.pop();
        if(grid->nodes[cur].status == TERM || grid->nodes[cur].status == PATH)
        {
            dst = cur;
            break;
        }
        for(auto& to : grid->adj[cur])
        {
            if(!used[to])
            {
                used[to] = true;
                queue.push(to);
                parents[to] = cur;
            }
        }
    }

    std::vector<int> path;
    grid->nodes[src].status = TERM;
    if(dst != -1)
    {
        for(int i = dst; i != -1; i = parents[i])
        {
            path.push_back(i);
            if(grid->nodes[i].status == NONE)
                grid->nodes[i].status = PATH;
        }
    }
    return path;
    
}

Graph* createSteinerTree(std::vector<Point>& points)
{
    std::sort(points.begin(), points.end());
    Point maxPoint = *(points.end() - 1);
    int size = std::max(maxPoint.x, maxPoint.y) + 1;

    
    Graph* grid = createGridGraph(size);
    auto ps = grid->getPointsId(points);
    

    Graph* steiner = new Graph();
    for(auto& p : ps)
    {
        auto path = addTerminal(grid, p);

        for(int i = 0; i < path.size(); i++)
        {
            int idx1 = steiner->addNode(grid->nodes[path[i]]);
            int idx2 = steiner->addNode(grid->nodes[path[i-1]]);
            if(i == 0)
                continue;
            steiner->addEdge(idx1, idx2);
        }
    } 

    delete grid;

    return steiner;
}


enum Layer
{
    PINS,
    M2,
    M3, 
    M2_M3,
    PINS_M2,
    PINS_M3
};

std::unordered_map<int, const char*> layers = {
    {PINS, "pins"},
    {M2, "m2"},
    {M3, "m3"},
    {M2_M3, "m2_m3"},
    {PINS_M2, "pins_m2"},
    {PINS_M3, "pins_m3"}};


struct OutPoint
{
    int x;
    int y;
    int layer;

    friend std::ostream& operator<<(std::ostream& os, const OutPoint& point);
};

std::ostream& operator<<(std::ostream& os, const OutPoint& point)
{
    os  << "<point x=\"" << point.x << "\" y=\"" << point.y 
        << "\" layer=\"" << layers[point.layer] << "\" />";
    return os;
}

struct Segment
{
    Point a;
    Point b;
    int layer;
    
    friend std::ostream& operator<<(std::ostream& os, const Segment& segment);
};

std::ostream& operator<<(std::ostream& os, const Segment& segment)
{
    os  << "<segment x1=\"" << segment.a.x << "\" y1=\"" << segment.a.y
        << "\" x2=\"" << segment.b.x << "\" y2=\"" << segment.b.y
        << "\" layer=\"" << layers[segment.layer] << "\" />";
    return os;
}


class OutputGenerator
{
public:
    void feed(Graph* g); 
    void generate(std::ifstream& input, std::ofstream& output);
private:
    std::vector<Segment> segments;
    std::vector<OutPoint> points;
};


void OutputGenerator::feed(Graph* g)
{
    int xsz = 0;
    int ysz = 0;

    for(auto& n : g->nodes)
    {
        if(xsz < n.x)
            xsz = n.x;
        if(ysz < n.y)
            ysz = n.y;
    }
    xsz++;
    ysz++;
    std::vector< std::vector<int> > vmap(ysz, std::vector<int> (xsz));
    std::vector< std::vector<int> > hmap(ysz, std::vector<int> (xsz));
  
    for(int node1 = 0; node1 < g->adj.size(); node1++)
    {
        for(auto& node2 : g->adj[node1])
        {
            if(std::abs(g->nodes[node1].x - g->nodes[node2].x) == 1)
            {
                hmap[g->nodes[node1].y][g->nodes[node1].x] = 1;
                hmap[g->nodes[node2].y][g->nodes[node2].x] = 1;
            }
            if(std::abs(g->nodes[node1].y - g->nodes[node2].y) == 1)
            {
                vmap[g->nodes[node1].y][g->nodes[node1].x] = 1;
                vmap[g->nodes[node2].y][g->nodes[node2].x] = 1;
            }
        }
    }
    
    for(int i = 0; i < ysz; i++)
    {
        for(int j = 0; j < xsz; j++)
        {
            int start = j;
            int end = start;
            if(hmap[i][j] == 1)
            {
                while(hmap[i][j++] == 1)
                {
                    end = j - 1;
                    if(j == xsz)
                        break;
                }

            }
            if(start != end)
                segments.push_back({ {start, i}, {end, i}, M2 });
        }
    }
     
    for(int j = 0; j < xsz; j++)
    {
        for(int i = 0; i < ysz; i++)
        {
            int start = i;
            int end = start;
            if(vmap[i][j] == 1)
            {
                while(vmap[i++][j] == 1)
                {
                    end = i - 1;
                    if(i == ysz)
                        break;
                }
                
            }
            if(start != end)
                segments.push_back({ {j, start}, {j, end}, M3 });
        }
    }
    
    std::vector<bool> external(g->size);

    for(auto& segment : segments)
    {
        int idx1 = g->find(segment.a);
        int idx2 = g->find(segment.b);
        external[idx1] = true;
        external[idx2] = true;
    }

    for(int i = 0; i < g->size; i++)
    {
        auto node = g->nodes[i];
        if(node.status == TERM)
            points.push_back({node.x, node.y, PINS_M2});
        if(hmap[node.y][node.x] == 1 && vmap[node.y][node.x] == 1)
            points.push_back({node.x, node.y, M2_M3});
        if(hmap[node.y][node.x] == 0 && vmap[node.y][node.x] == 1 && !external[i])
        {
            points.push_back({node.x, node.y, M2});
            points.push_back({node.x, node.y, M2_M3});
        }
    }

}

void OutputGenerator::generate(std::ifstream& input, std::ofstream& output)
{
    std::string tmp;
    while(std::getline(input, tmp))
    {
        if(tmp.find("</net") != std::string::npos)
        {
            for(auto& p : points)
                output << p << "\n";
            for(auto& s: segments)
                output << s << "\n";
        }
        output << tmp << "\n";
    }

}

int main(int argc, char* argv[])
{   
    if(argc == 1)
    {
        std::cerr << "No input file\n";
        std::exit(1);
    }

    std::string inName = std::string(argv[1]);
    std::ifstream input(inName);
    if(!input.is_open())
    {
        std::cerr << "Cannot open " << inName << "\n";
        std::exit(1);
    }
    std::string outName(inName);
    outName.replace(outName.begin() + outName.find(".xml"),
                    outName.end(), "_out.xml");

   
    auto points = parsePoints(input);
    Graph* steiner = createSteinerTree(points);
    
    input.clear();
    input.seekg(0);
    std::ofstream output(outName);  

    OutputGenerator generator;
    generator.feed(steiner);
    generator.generate(input, output);

    delete  steiner;
    return 0;
}


