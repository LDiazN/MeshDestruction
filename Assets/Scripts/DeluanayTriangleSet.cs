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

        public DeluanayTriangle(int v0, int v1, int v2, int adj0, int adj1, int adj2)
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

    private uint TriangleCount { get { return (uint) _triangles.Count / 3; } }

    public DeluanayTriangleSet(uint nPoints)
    {
        // Init containers
        _vertexPositions = new Vector2[nPoints + 3];
        _triangles = new List<int>();
        _adjacentTriangles = new List<int>();
        _recentlyAddedTriangles = new Stack<int>();
        InitFirstTriangle();
    }

    public void AddPoint(in Vector2 point)
    {
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

        _recentlyAddedTriangles.Push(t0Index);
        _recentlyAddedTriangles.Push(t1Index);
        _recentlyAddedTriangles.Push(t2Index);

        // TODO we still need to check deluanay condition
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
        _recentlyAddedTriangles.Push(0);
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

}
