using System.Collections;
using System.Collections.Generic;
using UnityEngine;

/// <summary>
/// Set of triangles data structure required to implement deluanay triangulation
/// </summary>
public class DeluanayTriangleSet
{
    /// <summary>
    /// Helper class to move data about triangles from here to there
    /// </summary>
    public struct DeluanayTriangle
    {
        public int[] vertices;
        public int[] adjacents;

        public DeluanayTriangle(int v0 = 0, int v1 = 0, int v2 = 0, int adj0 = -1, int adj1 = -1, int adj2 = -1)
        {
            vertices = new int[3] { v0, v1, v2 };
            adjacents = new int[3] { adj0, adj1, adj2 };
        }

    }

    /// <summary>
    /// Helper class to describe a triangle's edge
    /// </summary>
    public struct DeluanayTriangleEdge
    {
        int triangleIndex; // Which triangle this edge belongs to
        int edgeIndex; // Index of this edge inside the triangle 
        int vertexA, vertexB; // Indices of vertices the triangle vertex
    }

    public Vector2[] VertexPositions { get { return _vertexPositions; } }

    /// <summary>
    /// Store actual positions of mesh vertices
    /// </summary>
    private Vector2[] _vertexPositions;
    
    /// <summary>
    /// Array of triangle indices
    /// </summary>
    public List<int> Triangles { get { return _triangles; } }

    /// <summary>
    /// Triangles description, in the same format as a mesh. Three consecutive numbers
    /// enumerate the vertices of a triangle
    /// </summary>
    private List<int> _triangles;

    /// <summary>
    /// Holds the index of triangles adjacent to current triangle. 
    /// </summary>
    private List<int> _adjacentTriangles;

    /// <summary>
    /// Index of last added triangle
    /// </summary>
    private int _lastAddedTriangle;

    /// <summary>
    /// Current amount of points in vertexPositions array
    /// </summary>
    private uint _currentPointCount = 0; 

    /// <summary>
    /// Stack of recently added triangles, useful for deluanay edge swapping propagation
    /// </summary>
    private Stack<int> _recentlyAddedTriangles;

    /// <summary>
    /// How many triangles are currently stored
    /// </summary>
    private uint TriangleCount { get { return (uint) _triangles.Count / 3; } }

    /// <summary>
    /// Array of holes indices, we use this to control which points belong to which hole
    /// </summary>
    private List<int>[] _holesIndices; 

    public int HoleCount { get { return _holesIndices.Length; } }

    public DeluanayTriangleSet(uint nPoints, List<int> holeSizes)
    {
        // Init containers
        var nHoles = holeSizes.Count;

        // Sum size of each hole
        var sumOfHoleSizes = 0;

        foreach (var size in holeSizes)
            sumOfHoleSizes += size;

        _vertexPositions = new Vector2[nPoints + 3 + sumOfHoleSizes];
        _triangles = new List<int>();
        _adjacentTriangles = new List<int>();
        _recentlyAddedTriangles = new Stack<int>();
        InitFirstTriangle();

        // Create holes
        _holesIndices = new List<int>[nHoles];
        for (int i = 0; i < nHoles; i++)
            _holesIndices[i] = new List<int>();
    }

    /// <summary>
    /// Try to add a point to the triangulation. This function is guaranteed to keep a valid triangulation
    /// after each point is added
    /// </summary>
    /// <param name="point">new point to add</param>
    /// <param name="holeNum">if this point belongs to a hole. If < 0, it means that it isn't, and otherwise ti belongs to a hole with index `holeNum`</param>
    public void AddPoint(in Vector2 point, int holeNum = -1)
    {
        Debug.Assert(holeNum < HoleCount, "Hole number out of range");

        // Find a triangle containing this point, it must exists before we can proceed
        DeluanayTriangle? maybeTriangle;
        int t0Index;
        (maybeTriangle, t0Index) = FindContainingTriangle(point);
        DeluanayTriangle T0;
        // Sanity check
        if (maybeTriangle == null)
        {
            Debug.LogError($"Could not find a containing triangle for point: {point}");
            return;
        }
        else
            T0 = (DeluanayTriangle)maybeTriangle;

        // Add point 
        int newPointIndex = (int)_currentPointCount;
        _vertexPositions[_currentPointCount++] = point;

        // Add index to hole if necessary
        if (holeNum >= 0)
            _holesIndices[holeNum].Add(newPointIndex);

        // Create two triangles beside this
        int t1Index = (int)TriangleCount;
        int t2Index = t1Index + 1;

        // Create T1
        DeluanayTriangle T1 = new DeluanayTriangle(newPointIndex, T0.vertices[0], T0.vertices[1], t2Index, T0.adjacents[0], t0Index);

        // Create T2
        DeluanayTriangle T2 = new DeluanayTriangle(newPointIndex, T0.vertices[2], T0.vertices[0], t0Index, T0.adjacents[2], t1Index);

        // Update T0
        T0.vertices[0] = newPointIndex;
        T0.adjacents[0] = t1Index;
        T0.adjacents[2] = t2Index;

        // Add elements so it doesn't rise out of bounds
        _adjacentTriangles.Add(0);
        _adjacentTriangles.Add(0);
        _adjacentTriangles.Add(0);

        _triangles.Add(0);
        _triangles.Add(0);
        _triangles.Add(0);

        _adjacentTriangles.Add(0);
        _adjacentTriangles.Add(0);
        _adjacentTriangles.Add(0);

        _triangles.Add(0);
        _triangles.Add(0);
        _triangles.Add(0);

        // Set data of actual triangles
        SetDataOfTriangle(t1Index, T1);
        SetDataOfTriangle(t2Index, T2);
        SetDataOfTriangle(t0Index, T0);

        // now we need to update neighbors of t1 and t2 so they know who is their new neighbor
        if (T1.adjacents[1] != -1)
            for (int i = 0; i < 3; i++)
            {
                var adjI = 3 * T1.adjacents[1] + i;
                if (_adjacentTriangles[adjI] == t0Index)
                    _adjacentTriangles[adjI] = t1Index;
            }
        if (T2.adjacents[1] != -1)
            for (int i = 0; i < 3; i++)
            {
                var adjI = 3 * T2.adjacents[1] + i;
                if (_adjacentTriangles[adjI] == t0Index)
                    _adjacentTriangles[adjI] = t2Index;
            }

        if (T0.adjacents[1] != -1)
            _recentlyAddedTriangles.Push(t0Index);
        if (T1.adjacents[1] != -1)
            _recentlyAddedTriangles.Push(t1Index);
        if (T2.adjacents[1] != -1)
            _recentlyAddedTriangles.Push(t2Index);

        _lastAddedTriangle = t2Index;
        PropagateEdgeSwapping();
    }

    /// <summary>
    /// Init state of this class by adding the initial massive triangle
    /// </summary>
    private void InitFirstTriangle()
    {
        _vertexPositions[0] = new Vector2(-10, -10);
        _vertexPositions[1] = new Vector2(10, -10);
        _vertexPositions[2] = new Vector2(0, 10);
        _triangles.Add(0);
        _triangles.Add(2);
        _triangles.Add(1);
        _adjacentTriangles.Add(-1);
        _adjacentTriangles.Add(-1);
        _adjacentTriangles.Add(-1);
        _lastAddedTriangle = 0;
        _currentPointCount = 3;
    }
    
    /// <summary>
    /// Find a triangle containing the specified point if it exists.
    /// </summary>
    /// <param name="point">Point you want to check if it's contained by a triangle</param>
    /// <returns>A triangle if the specified point is contained inside some known triangle, null otherwise</returns>
    public (DeluanayTriangle?, int) FindContainingTriangle(in Vector2 point)
    {
        int nextTriangle = _lastAddedTriangle;

        DeluanayTriangle triangle = new DeluanayTriangle(0, 0, 0, -1, -1, -1);
        DeluanayTriangle? result = null;
        int resultIndex = -1;
        Vector3 init;
        Vector3 end;
        Vector3 p = new Vector3(point.x, point.y, 0);

        // find triangle containing this triangle set 
        while (nextTriangle != -1)
        {
            // Search edge with negative dot product
            GetDataOfTriangle(nextTriangle, ref triangle);
            for(int i = 0; i < 3; i++)
            {
                init = _vertexPositions[triangle.vertices[i]];
                end = _vertexPositions[triangle.vertices[(i+1) % 3]];
                var cross = Vector3.Cross(end - init, p - init);

                // Check if point is outside of this side.
                if (cross.z > 0)
                {
                    nextTriangle = triangle.adjacents[i];
                    break;
                }

                // If point inside all sides, then it's inside of this triangle
                if (i == 2)
                {
                    resultIndex = nextTriangle;
                    nextTriangle = -1;
                    result = triangle;
                }
            }
        }

        return (result, resultIndex);
    }

    /// <summary>
    /// Utility function to retrieve triangle data from internal arrays
    /// </summary>
    /// <param name="triangleIndex"> Index of triangle whose data you want to retrieve </param>
    /// <param name="triangle"> Where to store reulsts</param>
    void GetDataOfTriangle(int triangleIndex, ref DeluanayTriangle triangle)
    {
        if (triangleIndex > TriangleCount)
        {
            Debug.LogError($"Index: {triangleIndex} is out of range of the current triangle count: {TriangleCount}");
            return;
        }

        triangle.vertices[0] = _triangles[triangleIndex * 3];
        triangle.vertices[1] = _triangles[triangleIndex * 3 + 1];
        triangle.vertices[2] = _triangles[triangleIndex * 3 + 2];

        triangle.adjacents[0] = _adjacentTriangles[triangleIndex * 3];
        triangle.adjacents[1] = _adjacentTriangles[triangleIndex * 3 + 1];
        triangle.adjacents[2] = _adjacentTriangles[triangleIndex * 3 + 2];
    }

    /// <summary>
    /// Utility function to just update the values of a triangle of index `triangleIndex` 
    /// </summary>
    /// <param name="triangleIndex">index of triangle to update</param>
    /// <param name="triangle">new data to update</param>
    void SetDataOfTriangle(int triangleIndex, in DeluanayTriangle triangle)
    {
        // Sanity check
        if (triangleIndex > TriangleCount)
        {
            Debug.LogError($"Index: {triangleIndex} ({3 * triangleIndex}) is out of range of the current triangle count: {TriangleCount}");
            return;
        }

        _triangles[triangleIndex * 3] = triangle.vertices[0];
        _triangles[triangleIndex * 3 + 1] = triangle.vertices[1];
        _triangles[triangleIndex * 3 + 2] = triangle.vertices[2];

        _adjacentTriangles[triangleIndex * 3] = triangle.adjacents[0];
        _adjacentTriangles[triangleIndex * 3 + 1] = triangle.adjacents[1];
        _adjacentTriangles[triangleIndex * 3 + 2] = triangle.adjacents[2];
    }

    /// <summary>
    /// Check deluanay condition in the triangle formed by v0,v1,v2 and the point p. This is, check if p is outside the 
    /// circumcircle formed by the triangle, specified in CCW order. Computed using the formula from here: 
    /// https://en.wikipedia.org/wiki/Delaunay_triangulation
    /// </summary>
    /// <param name="v0">First Vertex</param>
    /// <param name="v1">Second Vertex</param>
    /// <param name="v2">Third Vertex</param>
    /// <param name="p">Point to check if it's inside the triangle</param>
    /// <returns>true triangle and point match deluanay condition. This is, when p is outside the circumcircle. false otherwise</returns>
    static bool CheckDeluanayCondition(in Vector2 v0, in Vector2 v1, in Vector2 v2, in Vector2 p)
    {
        var A = v0;
        var B = v1;
        var C = v2;
        var D = p;

        var mat = new Matrix4x4(
                new Vector4(A.x, B.x, C.x, D.x), 
                new Vector4(A.y, B.y, C.y, D.y),
                new Vector4(A.x*A.x + A.y*A.y, B.x*B.x + B.y*B.y, C.x*C.x + C.y*C.y, D.x*D.x + D.y*D.y),
                Vector4.one
            );

        return mat.determinant > 0;
    }

    /// <summary>
    /// Try to propagate edge swapping according to the stack of recently added triangles
    /// </summary>
    void PropagateEdgeSwapping()
    {
        DeluanayTriangle triangle = new DeluanayTriangle(0,0,0, -1, -1, -1);
        DeluanayTriangle oppositeTriangle = new DeluanayTriangle(0,0,0, -1, -1, -1);

        while (_recentlyAddedTriangles.Count != 0)
        {
            var nextTriangle = _recentlyAddedTriangles.Pop(); // Ti
            GetDataOfTriangle(nextTriangle, ref triangle);
            var Ta = triangle.adjacents[1];
            if (Ta == -1)
                continue; // there's no opposite triangle, continue
            GetDataOfTriangle(Ta, ref oppositeTriangle); // Ta

            // Try to check if this triangle is 
            if (CheckDeluanayCondition(
                _vertexPositions[oppositeTriangle.vertices[0]], 
                _vertexPositions[oppositeTriangle.vertices[1]], 
                _vertexPositions[oppositeTriangle.vertices[2]], 
                _vertexPositions[triangle.vertices[0]])
            )
                continue;

            // Now that we know that we need an edge swapping, we have to gather other edges of the 
            // adjacent triangles to push them in the stack and propagate edge swapping
            int commonEdgeIndex = -1; // Ea
            for (int i = 0; i < 3; i++)
                if (oppositeTriangle.adjacents[i] == nextTriangle)
                    commonEdgeIndex = i;

            int Va = (commonEdgeIndex + 2) % 3;
            int Vi = 0;

            // Sanity check
            Debug.Assert(commonEdgeIndex != -1, "Inconsistent adjacent triangle, one of the sides should be nextTriangle");

            // Stack neighboring triangles
            var nextTriangleToPush = oppositeTriangle.adjacents[(commonEdgeIndex + 1) % 3];
            if (nextTriangleToPush != -1)
                _recentlyAddedTriangles.Push(nextTriangleToPush);
            nextTriangleToPush = oppositeTriangle.adjacents[(commonEdgeIndex + 2) % 3];
            if (nextTriangleToPush != -1)
                _recentlyAddedTriangles.Push(nextTriangleToPush);

            // Perform edge swapping
            triangle.vertices[(Vi + 1) % 3] = oppositeTriangle.vertices[Va];
            oppositeTriangle.vertices[commonEdgeIndex] = triangle.vertices[Vi];
            oppositeTriangle.adjacents[commonEdgeIndex] = triangle.adjacents[Vi];
            triangle.adjacents[Vi] = Ta;
            triangle.adjacents[(Vi + 1) % 3] = oppositeTriangle.adjacents[Va];
            oppositeTriangle.adjacents[Va] = nextTriangle;


            SetDataOfTriangle(Ta, oppositeTriangle);
            SetDataOfTriangle(nextTriangle, triangle);

            // Swap edges with neighbor
            var Tb = triangle.adjacents[(Vi + 1) % 3];
            var Tc = oppositeTriangle.adjacents[commonEdgeIndex];

            // Search neighbor with out dated adjacent and update it
            if (Tb != -1)
                for (int i = 0; i < 3; i++)
                {
                    var index = 3 * Tb + i;
                    if (_adjacentTriangles[index] == Ta)
                        _adjacentTriangles[index] = nextTriangle;
                }
            // Search neighbor with out dated adjacent and update it
            if (Tc != -1)
                for (int i = 0; i < 3; i++)
                {
                    var index = 3 * Tc + i;
                    if (_adjacentTriangles[index] == nextTriangle)
                        _adjacentTriangles[index] = Ta;
                }
        }
    }

    public void FixConstrainedEdges()
    {
        DeluanayTriangle triangle = new DeluanayTriangle(0,0,0, -1, -1, -1);
        bool[] intersectedEdges = new bool[3] {false, false, false};

        foreach(var hole in _holesIndices)
        {
            for(int holeIndex = 0; holeIndex < hole.Count; holeIndex += 1)
            {
                var constEdgeStart = hole[holeIndex];
                var constEdgeEnd = hole[(holeIndex+1) % hole.Count];

                var constEdgeStartPoint = _vertexPositions[constEdgeStart];
                var constEdgeEndPoint = _vertexPositions[constEdgeEnd];

                // Search starting triangle. This is a triangle such that 
                // one of the vertices is the current edgeStart, and other is the edgeEnd
                // or one of the vertices is edgeStart and one of the sides intersects this edge

                // We want the first tri to iterate over in this variable
                var nextTri = -1;
                List<(int, int)> intersectingEdges = new List<(int, int)>();
                for(int tri = 0; tri < TriangleCount; tri ++)
                {
                    GetDataOfTriangle(tri, ref triangle);

                    var vi = -1;

                    // Skip if none of them is edgeStart
                    if (triangle.vertices[0] == constEdgeStart)
                        vi = 0;
                    else if (triangle.vertices[1] == constEdgeStart)
                        vi = 1;
                    else if (triangle.vertices[2] == constEdgeStart)
                        vi = 2;
                    else
                        continue;

                    // Sanity check
                    Debug.Assert(vi != -1, "Should only get here if you found vi");

                    var vi1 = (vi + 1) % 3;
                    var vi2 = (vi + 2) % 3;

                    var viVert = _vertexPositions[triangle.vertices[vi]];
                    var vi1Vert = _vertexPositions[triangle.vertices[vi1]];
                    var vi2Vert = _vertexPositions[triangle.vertices[vi2]];
                    var edgeEndVert = _vertexPositions[constEdgeEnd];

                    var E1 = vi1Vert - viVert;
                    var E2 = viVert - vi2Vert;

                    // If in the left side of E1, E2 at the same time...
                    // TODO double check this innequalities
                    if (Cross2D(E1, edgeEndVert - viVert) <= 0 && Cross2D(E2, edgeEndVert - vi2Vert) <= 0)
                    {
                        nextTri = tri;
                        break;
                    }
                }

                // Sanity check
                Debug.Assert(nextTri != -1, "There should be at the list one triangle intersecting");

                // While next tri is not the one containing edgeEnd
                GetDataOfTriangle(nextTri, ref triangle);
                while (triangle.vertices[0] != constEdgeEnd && triangle.vertices[1] != constEdgeEnd && triangle.vertices[2] != constEdgeEnd)
                {
                    intersectedEdges[0] = intersectedEdges[1] = intersectedEdges[2] = false;

                    var Et = -1;
                    for(int Ei = 0; Ei < 3; Ei++)
                    {
                        // Check if this edge is in right side (opposite side as before)
                        var edgeStartVert = _vertexPositions[triangle.vertices[Ei]];
                        var edgeEndVert = _vertexPositions[triangle.vertices[(Ei+1) % 3]];

                        // Remember that edgeEndPoint is related to the hole, edgeEndVert is related
                        // to the current edge in this triangle
                        if (Cross2D(edgeEndVert - edgeStartVert, constEdgeEndPoint - edgeStartVert) <= 0)
                            Et = Ei;
                        else
                            continue;

                        // Check intersection between this edge and the constrained edge
                        if (!CheckSegmentIntersection(edgeStartVert, edgeEndVert, constEdgeStartPoint, constEdgeEndPoint))
                            continue;

                        // Since they do intersect, we register edges, set next tri
                        intersectingEdges.Add((triangle.vertices[holeIndex], triangle.vertices[(Ei+1) % 3]));
                        intersectedEdges[Ei] = true;
                        intersectedEdges[(Ei + 1) % 3] = true;
                        nextTri = triangle.adjacents[Ei];
                        break;
                    }

                    // If no intersection...
                    if (!intersectedEdges[0] && !intersectedEdges[1] && !intersectedEdges[2])
                        nextTri = triangle.adjacents[Et];


                    GetDataOfTriangle(nextTri, ref triangle);
                }
            }
        }
    }

    /// <summary>
    /// Check if the two line segments specified by starting point p1, end point p2,
    /// starting point q1, and end point q2 intersect in a single point. Implemented according 
    /// to this post: https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
    /// </summary>
    /// <param name="p1">Starting point of line segment 1</param>
    /// <param name="p2">End point of line segment 1</param>
    /// <param name="q1">Starting point of line segment 2</param>
    /// <param name="q2">End point of line segment 2</param>
    /// <returns></returns>
    private static bool CheckSegmentIntersection(in Vector2 p1, in Vector2 p2, in Vector2 q1, in Vector2 q2)
    {

        var p = p1;
        var r = p2 - p1;

        var q = q1;
        var s = q2 - q1;

        var rXs = Cross2D(r, s);
        var q_pXr = Cross2D(q - p, r);

        const float delta = 0.0001f;


        if (
            Mathf.Abs(rXs) < delta && Mathf.Abs(q_pXr) < delta || // Collinear
            Mathf.Abs(rXs) < delta  // Parallel

            ) 
            return false;

        var t = Cross2D(q - p, (1.0f / (Cross2D(r, s)) * s));
        var u = Cross2D(q - p, (1/0f / Cross2D(r,s)) * r);

        // If values are in valid range, they intersect. Otherwise, 
        // they are not parallel but they don't intersect
        return 0 <= t && t <= 1 && 0 <= u && u <= 1;
    }

    private static float Cross2D(in Vector2 a, in Vector2 b) => a.x * b.y - a.y * b.x;
}
