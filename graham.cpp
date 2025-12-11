#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <algorithm>

using namespace std;
using namespace std::chrono;

const double EPSILON = 1e-9;

// Metrics structure
struct Metrics {
    double preprocessTime = 0;
    double sortingTime = 0; 
    double algorithmTime = 0;
    double totalRuntime = 0;
    int orientationTests = 0;
    int comparisons = 0;
    int inputSizeAfterDuplicates = 0;
    int hullSize = 0;
    double peakMemoryKB = 0;
    int maxStackSize = 0;
};

Metrics metrics;

// Point structure
struct Point {
    double x, y;
    Point(double x = 0, double y = 0) : x(x), y(y) {}
};

// File Validation Class
class FileValidation {
public:
    bool checkFileExists(const string& filename) {
        ifstream file(filename);
        if (!file.is_open()) {
            return false;
        }
        file.close();
        return true;
    }
    
    bool checkFileReadable(const string& filename) {
        ifstream file(filename);
        if (!file.good()) {
            return false;
        }
        string line;
        try {
            if (!getline(file, line)){
                file.close();
                return false;
            }
        } catch (...) {
            file.close();
            return false;
        }
        if (file.bad() || file.fail()) {
            file.close();
            return false;
        }
        file.close();
        return true;
    }
    
    bool checkFileEmpty(const string& filename) {
        ifstream file(filename);
        string content, line;
        while (getline(file, line)) {
            for (char c : line) {
                if (!isspace(c)) {
                    file.close();
                    return true;
                }
            }
        }
        file.close();
        return false;
    }
    
    bool validateCoordinateFormat(const string& line) {
        if (line.empty()) return false;
        
        if (line.find(',') == string::npos) return false;
        
        stringstream ss(line);
        string token;
        vector<string> tokens;
        
        while (getline(ss, token, ',')) {
            token.erase(0, token.find_first_not_of(" \t\r\n"));
            token.erase(token.find_last_not_of(" \t\r\n") + 1);
            if (!token.empty()) {
                tokens.push_back(token);
            }
        }
        
        if (tokens.empty()) return false;
        
        for (const string& t : tokens) {
            bool hasDigit = false;
            bool hasDecimal = false;
            
            for (size_t i = 0; i < t.length(); i++) {
                char c = t[i];
                if (isdigit(c)) {
                    hasDigit = true;
                } else if (c == '.') {
                    if (hasDecimal) return false;
                    hasDecimal = true;
                } else if (c == '-' || c == '+') {
                    if (i != 0) return false;
                } else if (c == 'e' || c == 'E') {
                    if (!hasDigit) return false;
                    hasDigit = false;
                } else {
                    return false;
                }
            }
            if (!hasDigit) return false;
        }
        return true;
    }
    
    bool checkDimensionsCorrect(const vector<string>& lines) {
        for (const string& line : lines) {
            if (line.empty()) continue;
            
            stringstream ss(line);
            string token;
            int count = 0;
            
            while (getline(ss, token, ',')) {
                token.erase(0, token.find_first_not_of(" \t\r\n"));
                token.erase(token.find_last_not_of(" \t\r\n") + 1);
                if (!token.empty()) count++;
            }
            
            if (count != 2) {
                return false;
            }
        }
        return true;
    }
    
    int parseFile(const string& filename, vector<Point>& points) {
        ifstream file(filename);
        string line;
        int count = 0;
        
        while (getline(file, line)) {
            if (line.empty()) continue;
            
            stringstream ss(line);
            string xStr, yStr;
            
            if (getline(ss, xStr, ',') && getline(ss, yStr, ',')) {
                xStr.erase(0, xStr.find_first_not_of(" \t\r\n"));
                xStr.erase(xStr.find_last_not_of(" \t\r\n") + 1);
                yStr.erase(0, yStr.find_first_not_of(" \t\r\n"));
                yStr.erase(yStr.find_last_not_of(" \t\r\n") + 1);
                
                try {
                    double x = stod(xStr);
                    double y = stod(yStr);
                    
                    if (isnan(x) || isnan(y) || isinf(x) || isinf(y)) {
                        continue;
                    }
                    
                    points.push_back(Point(x, y));
                    count++;
                } catch (...) {
                    continue;
                }
            }
        }
        file.close();
        return count;
    }
};

// Duplicate Filtering Class
class DuplicateFiltering {
public:
    bool removeDuplicates(vector<Point>& points, int& numPoints) {
        if (points.empty()) return false;
        
        vector<Point> unique;
        int counter = 0;
        
        for (int i = 0; i < numPoints; i++) {
            bool isDuplicate = false;
            
            for (int j = 0; j < counter; j++) {
                if (fabs(points[i].x - unique[j].x) < EPSILON && 
                    fabs(points[i].y - unique[j].y) < EPSILON) {
                    isDuplicate = true;
                    break;
                }
            }
            
            if (!isDuplicate) {
                unique.push_back(points[i]);
                counter++;
            }
        }
        
        points = unique;
        numPoints = counter;
        metrics.inputSizeAfterDuplicates = counter;
        return counter >= 3;
    }
    
    bool checkCollinear(const vector<Point>& points, int numPoints) {
        if (numPoints < 3) return true;
        
        bool isCollinear = true;
        
        double dx = points[1].x - points[0].x;
        double dy = points[1].y - points[0].y;
        double slope1;
        
        if (fabs(dx) < EPSILON) {
            slope1 = INFINITY;
        } else {
            slope1 = fabs(dy / dx);
        }
        
        for (int i = 2; i < numPoints; i++) {
            double dx2 = points[i].x - points[0].x;
            double dy2 = points[i].y - points[0].y;
            double slope2;
            
            if (fabs(dx2) < EPSILON) {
                slope2 = INFINITY;
            } else {
                slope2 = fabs(dy2 / dx2);
            }
            
            if (fabs(slope1 - slope2) > EPSILON && 
                !(isinf(slope1) && isinf(slope2))) {
                isCollinear = false;
                break;
            }
        }
        
        return isCollinear;
    }
};

// Stack Class
class Stack {
private:
    vector<Point> arr;
    int top;
    
public:
    Stack(int capacity) {
        arr.resize(capacity);
        top = -1;
    }
    
    void push(const Point& p) {
        if (top + 1 >= (int)arr.size()) {
            arr.resize(arr.size() * 2);
        }
        arr[++top] = p;
        
        // Track maximum stack size
        if (top + 1 > metrics.maxStackSize) {
            metrics.maxStackSize = top + 1;
        }
    }
    
    Point pop() {
        if (top < 0) {
            return Point(0, 0);
        }
        return arr[top--];
    }
    
    Point peek() {
        if (top < 0) {
            return Point(0, 0);
        }
        return arr[top];
    }
    
    Point nextToTop() {
        if (top < 1) {
            return Point(0, 0);
        }
        return arr[top - 1];
    }
    
    int size() {
        return top + 1;
    }
    
    vector<Point> getAll() {
        vector<Point> result;
        for (int i = 0; i <= top; i++) {
            result.push_back(arr[i]);
        }
        return result;
    }
};

// Graham Scan Algorithm Class
class GrahamScan {
private:
    Point anchor;
    
public:
    int findLowest(const vector<Point>& points, int numPoints) {
        int lowest = 0;
        
        for (int i = 1; i < numPoints; i++) {
            metrics.comparisons++;
            if (points[i].y < points[lowest].y - EPSILON) {
                lowest = i;
            } else if (fabs(points[i].y - points[lowest].y) < EPSILON) {
                metrics.comparisons++;
                if (points[i].x < points[lowest].x - EPSILON) {
                    lowest = i;
                }
            }
        }
        return lowest;
    }
    
    double polarAngle(const Point& target) {
        double dx = target.x - anchor.x;
        double dy = target.y - anchor.y;
        return atan2(dy, dx);
    }
    
    double euclideanDistance(const Point& p1, const Point& p2) {
        double dx = p2.x - p1.x;
        double dy = p2.y - p1.y;
        return sqrt(dx * dx + dy * dy);
    }
    
    void swapPoints(vector<Point>& points, int i, int j) {
        Point temp = points[i];
        points[i] = points[j];
        points[j] = temp;
    }
    
    void sortByPolarAngle(vector<Point>& points, int numPoints, int anchorIndex) {
        swapPoints(points, 0, anchorIndex);
        anchor = points[0];
        
        // Sorting with comparison counting
        for (int i = 1; i < numPoints - 1; i++) {
            for (int j = i + 1; j < numPoints; j++) {
                metrics.comparisons++;
                double angle1 = polarAngle(points[i]);
                double angle2 = polarAngle(points[j]);
                
                if (angle1 > angle2 + EPSILON) {
                    swapPoints(points, i, j);
                } else if (fabs(angle1 - angle2) < EPSILON) {
                    metrics.comparisons++;
                    double dist1 = euclideanDistance(anchor, points[i]);
                    double dist2 = euclideanDistance(anchor, points[j]);
                    if (dist1 > dist2 + EPSILON) {
                        swapPoints(points, i, j);
                    }
                }
            }
        }
    }
    
    int orientation(const Point& p, const Point& q, const Point& r) {
        metrics.orientationTests++;
        double val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
        
        if (fabs(val) < EPSILON) return 0;
        return (val > 0) ? 1 : 2;
    }
    
    void handleCollinearEndPoints(vector<Point>& points, int numPoints) {
        int i = numPoints - 2;
        
        while (i >= 0) {
            int orient = orientation(points[0], points[i], points[numPoints - 1]);
            if (orient == 0) {
                i--;
            } else {
                break;
            }
        }
        
        i++;
        
        if (i < numPoints - 1) {
            int left = i;
            int right = numPoints - 1;
            while (left < right) {
                swapPoints(points, left, right);
                left++;
                right--;
            }
        }
    }
    
    int computeConvexHull(vector<Point>& points, int numPoints, vector<Point>& hull) {
        if (numPoints < 3) return 0;
        
        int anchorIndex = findLowest(points, numPoints);

        auto sortStart = high_resolution_clock::now();
        sortByPolarAngle(points, numPoints, anchorIndex);
        handleCollinearEndPoints(points, numPoints);
        auto sortEnd = high_resolution_clock::now();
        metrics.sortingTime = duration_cast<microseconds>(sortEnd - sortStart).count() / 1000.0;
        
        Stack stack(numPoints);
        stack.push(points[0]);
        stack.push(points[1]);
        stack.push(points[2]);
        
        for (int i = 3; i < numPoints; i++) {
            metrics.comparisons++;
            while (stack.size() > 1) {
                Point top = stack.peek();
                Point nextToTop = stack.nextToTop();
                
                int orient = orientation(nextToTop, top, points[i]);
                if (orient != 2) {
                    stack.pop();
                } else {
                    break;
                }
            }
            stack.push(points[i]);
        }
        
        hull = stack.getAll();
        metrics.hullSize = hull.size();
        return hull.size();
    }
    
    void printResults(const vector<Point>& hull, int hullSize) {
        cout << "\nConvex Hull Vertices (in order):\n";
        for (int i = 0; i < hullSize; i++) {
            cout << "Point " << i << ": (" << fixed << setprecision(6) 
                 << hull[i].x << ", " << hull[i].y << ")\n";
        }
        cout << "Dimension of the convex hull: " << hullSize << "\n";
    }
    
    void writeResults(const vector<Point>& hull, int hullSize, const string& filename) {
        ofstream file(filename);
        if (!file.is_open()) {
            cout << "Error: Cannot open output file\n";
            return;
        }
        
        file << "Convex Hull Vertices (in order):\n";
        for (int i = 0; i < hullSize; i++) {
            file << "Point " << i << ": (" << fixed << setprecision(6) 
                 << hull[i].x << ", " << hull[i].y << ")\n";
        }
        file << "Dimension of the convex hull: " << hullSize << "\n";
        file.close();
        
        cout << "Results written to " << filename << "\n";
    }
};

// Ray Casting Algorithm for correctness verification
class CorrectnessVerification {
public:
    bool isPointInPolygon(const Point& p, const vector<Point>& hull) {
        int n = hull.size();
        int count = 0;
        
        for (int i = 0; i < n; i++) {
            int j = (i + 1) % n;
            
            if ((hull[i].y <= p.y && p.y < hull[j].y) || 
                (hull[j].y <= p.y && p.y < hull[i].y)) {
                
                double xIntersect = hull[i].x + (p.y - hull[i].y) * 
                                    (hull[j].x - hull[i].x) / (hull[j].y - hull[i].y);
                
                if (p.x < xIntersect) {
                    count++;
                }
            }
        }
        
        return (count % 2) == 1;
    }
    
    double verifyHull(const vector<Point>& allPoints, const vector<Point>& hull) {
        int outliers = 0;
        int total = allPoints.size();
        
        for (const Point& p : allPoints) {
            bool onHull = false;
            for (const Point& h : hull) {
                if (fabs(p.x - h.x) < EPSILON && fabs(p.y - h.y) < EPSILON) {
                    onHull = true;
                    break;
                }
            }
            
            if (!onHull && !isPointInPolygon(p, hull)) {
                outliers++;
            }
        }
        
        return 100.0 * (total - outliers) / total;
    }
};


void printMetrics(const string& filename) {
    ofstream file(filename, ios::app);
    cout << "\n========================================\n";
    cout << "PERFORMANCE METRICS\n";
    cout << "========================================\n";
    
    cout << "\nTime-Based Metrics:\n";
    cout << "  Preprocessing Time: " << fixed << setprecision(3) 
         << metrics.preprocessTime << " ms\n";
    cout << "  Sorting Time: " << metrics.sortingTime << " ms\n";  // NEW LINE
    cout << "  Algorithm Time (hull construction): " << metrics.algorithmTime << " ms\n";
    cout << "  Total Runtime: " << metrics.totalRuntime << " ms\n";
    
    cout << "\nOperation Count Metrics:\n";
    cout << "  Orientation Tests: " << metrics.orientationTests << "\n";
    cout << "  Number of Comparisons: " << metrics.comparisons << "\n";
    
    cout << "\nSpace Metrics:\n";
    cout << "  Peak Memory Usage: " << fixed << setprecision(2) 
         << metrics.peakMemoryKB << " KB\n";
    cout << "  Stack Maximum Size: " << metrics.maxStackSize << "\n";  // NEW LINE
    
    cout << "\nInput/Output Characterization:\n";
    cout << "  Input Size (after duplicate removal): " << metrics.inputSizeAfterDuplicates << "\n";
    cout << "  Hull Size: " << metrics.hullSize << "\n";
    
    cout << "========================================\n";

    // Write to file as well
    file << "\n========================================\n";
    file << "PERFORMANCE METRICS\n";
    file << "========================================\n";
    
    file << "\nTime-Based Metrics:\n";
    file << "  Preprocessing Time: " << fixed << setprecision(3) 
         << metrics.preprocessTime << " ms\n";
    file << "  Sorting Time: " << metrics.sortingTime << " ms\n";  // NEW LINE
    file << "  Algorithm Time (hull construction): " << metrics.algorithmTime << " ms\n";
    file << "  Total Runtime: " << metrics.totalRuntime << " ms\n";
    
    file << "\nOperation Count Metrics:\n";
    file << "  Orientation Tests: " << metrics.orientationTests << "\n";
    file << "  Number of Comparisons: " << metrics.comparisons << "\n";
    
    file << "\nSpace Metrics:\n";
    file << "  Peak Memory Usage: " << fixed << setprecision(2) 
         << metrics.peakMemoryKB << " KB\n";
    file << "  Stack Maximum Size: " << metrics.maxStackSize << "\n";  // NEW LINE
    
    file << "\nInput/Output Characterization:\n";
    file << "  Input Size (after duplicate removal): " << metrics.inputSizeAfterDuplicates << "\n";
    file << "  Hull Size: " << metrics.hullSize << "\n";
    
    file << "========================================\n";

    file.close();
}

void writeError(const string& filename, const string& error){
    ofstream file(filename);
    if (!file.is_open()) {
        cout << "Error: Cannot open output file\n";
        return;
    }
    file << error << "\n";
    file.close();
}

int main(int argc, char* argv[]) {
    auto totalStart = high_resolution_clock::now();
    
    string inputFilename, outputFilename;
    
    cout << "========================================\n";
    cout << "GRAHAM SCAN CONVEX HULL ALGORITHM\n";
    cout << "========================================\n";
    if (argc == 1){
        cout << "Enter input CSV filename: ";
        cin >> inputFilename;
        cout << "Enter output filename: ";
        cin >> outputFilename;
    }
    else if (argc == 3){
        inputFilename = argv[1];
        outputFilename = argv[2];
    }
    else{
        cerr << "Usage: Either just " << argv[0] << " or " << argv[0] << " (inputFilename) (outputFilename)\n";
        return 1;
    }
    
    FileValidation validator;
    
    if (!validator.checkFileExists(inputFilename)) {
        cout << "Error: File Not Found\n";
        writeError(outputFilename, "Error: File Not Found\n");
        return 1;
    }
    
    if (!validator.checkFileReadable(inputFilename)) {
        cout << "Error: Input/Output Error - Cannot read file\n";
        writeError(outputFilename, "Error: Input/Output Error - Cannot read file\n");
        return 1;
    }
    
    if (!validator.checkFileEmpty(inputFilename)) {
        cout << "Error: Empty File\n";
        writeError(outputFilename, "Error: Empty File\n");
        return 1;
    }
    
    ifstream file(inputFilename);
    vector<string> lines;
    string line;
    while (getline(file, line)) {
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        if (!line.empty()) lines.push_back(line);
    }
    file.close();
    
    if (lines.empty()) {
        cout << "Error: Empty File\n";
        writeError(outputFilename, "Error: Empty File\n");
        return 1;
    }
    
    if (!validator.checkDimensionsCorrect(lines)) {
        cout << "Error: Dimension Mismatch - Either Inconsistent Dimensions or non-2D Dimensions\n";
        writeError(outputFilename, "Error: Dimension Mismatch - Either Inconsistent Dimensions or non-2D Dimensions\n");
        return 1;
    }
    
    for (const string& l : lines) {
        if (!validator.validateCoordinateFormat(l)) {
            cout << "Error: Invalid Data Format - Non-numeric or unexpected coordinates\n";
            writeError(outputFilename, "Error: Invalid Data Format - Non-numeric or unexpected coordinates\n");
            return 1;
        }
    }
    vector<Point> points;
    int numPoints = validator.parseFile(inputFilename, points);
    
    if (numPoints == 0) {
        cout << "Error: No valid points found in file\n";
        writeError(outputFilename, "Error: No valid points found in file\n");
        return 1;
    }
    
    // Start preprocessing timer
    auto preprocessStart = high_resolution_clock::now();
    
    DuplicateFiltering filter;
    if (!filter.removeDuplicates(points, numPoints)) {
        cout << "Error: Insufficient Points - At least 3 unique points required\n";
        writeError(outputFilename, "Error: Insufficient Points - At least 3 unique points required\n");
        return 1;
    }
    
    if (filter.checkCollinear(points, numPoints)) {
        cout << "Error: Collinear Points - All points lie on a single line\n";
        writeError(outputFilename, "Error: Collinear Points - All points lie on a single line\n");
        return 1;
    }
    
    
    auto preprocessEnd = high_resolution_clock::now();
    metrics.preprocessTime = duration_cast<microseconds>(preprocessEnd - preprocessStart).count() / 1000.0;
    
    // Start algorithm timer
    auto algoStart = high_resolution_clock::now();
    
    GrahamScan graham;
    vector<Point> hull;
    int hullSize = graham.computeConvexHull(points, numPoints, hull);
    
    auto algoEnd = high_resolution_clock::now();
    metrics.algorithmTime = duration_cast<microseconds>(algoEnd - algoStart).count() / 1000.0;

    // Subtract the sorting time to get the actual algorithm time
    metrics.algorithmTime -= metrics.sortingTime;
    
    if (hullSize == 0) {
        cout << "Error: Failed to compute convex hull\n";
        writeError(outputFilename, "Error: Failed to compute convex hull\n");
        return 1;
    }
    
    // Calculate peak memory usage (points + hull + stack)
    metrics.peakMemoryKB = (sizeof(Point) * points.size() + 
                            sizeof(Point) * hull.size() + 
                            sizeof(Point) * metrics.maxStackSize) / 1024.0;
    
    // Correctness verification using Ray Casting
    CorrectnessVerification verifier;
    double accuracy = verifier.verifyHull(points, hull);
    
    graham.printResults(hull, hullSize);
    graham.writeResults(hull, hullSize, outputFilename);
    
    auto totalEnd = high_resolution_clock::now();
    metrics.totalRuntime = duration_cast<microseconds>(totalEnd - totalStart).count() / 1000.0;
    
    printMetrics(outputFilename);
    
    cout << "\nCorrectness Metrics:\n";
    cout << "  Ray Casting Verification Accuracy: " << fixed << setprecision(2) 
         << accuracy << "%\n";
    
    cout << "\nConvex hull computation completed successfully\n";
    
    return 0;
}