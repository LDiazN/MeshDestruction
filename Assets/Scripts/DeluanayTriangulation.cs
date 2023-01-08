using System.Collections;
using System.Collections.Generic;
using UnityEngine;

/// <summary>
/// Implements a constrained deluanay triangulation, based on this implementation
/// https://forum.unity.com/threads/programming-tools-constrained-delaunay-triangulation.1066148/
/// </summary>
public class DeluanayTriangulation 
{
    struct PointNormalizationResult
    {
        public Vector2[] normalizedPoints;
        public float height;
        public float width;
        public float maxX, minX;
        public float maxY, minY;
    }

    public struct TriangulationResult
    {
        public int[] triangles;
        public Vector2[] points;
    }

    private DeluanayTriangulation()
    {

    }

    /// <summary>
    /// This triangulation algorithm will leverage the algorithm in the following paper: https://cg.cs.uni-bonn.de/backend/v1/files/publications/klein-1996-construction.pdf
    /// </summary>
    /// <param name="pointsToTriangulate"></param>
    /// <param name="holes"></param>
    /// <param name="isCounterClockWise"></param>
    /// <returns>Triangulation result object specifying desired data</returns>
    public static TriangulationResult Triangulate(in Vector2[] pointsToTriangulate, in List<Vector2[]> holes, bool isCounterClockWise = true)
    {
        TriangulationResult result;
        // Create linked list with starting polygon
        LinkedList<int> polygon = new LinkedList<int>();

        if (isCounterClockWise)
            for (int i = 0; i < pointsToTriangulate.Length; i++)
                polygon.AddLast(i);
        else
            for (int i = pointsToTriangulate.Length-1; i > 0 ; i--)
                polygon.AddLast(i);

        List<Vector3Int> results = new List<Vector3Int>();
        TriangulatePSLG(polygon, pointsToTriangulate, ref results);

        // Build solution
        int[] triangles = new int[results.Count * 3];
        for (int i = 0; i < results.Count; i++)
        {
            triangles[i * 3] = results[i].x;
            triangles[i * 3 + 1] = results[i].y;
            triangles[i * 3 + 2] = results[i].z;
        }
        result.triangles = triangles;
        result.points = pointsToTriangulate;
        
        return result;
    }

    public static TriangulationResult Triangulate_old(in Vector2[] pointsToTriangulate, in List<Vector2[]> holes)
    {
        var pointNormalization = NormalizePoints(pointsToTriangulate);
        var pointBinGrid = new PointBinGrid(pointNormalization.normalizedPoints);
        
        // Count how many vertex per hole
        List<int> holeSizes = new();
        foreach (var hole in holes)
            holeSizes.Add(hole.Length);

        var triangleSet = new DeluanayTriangleSet((uint) pointsToTriangulate.Length, holeSizes); 

        // Normalize holes using the same normalization process for the input points
        var D = Mathf.Max(pointNormalization.width, pointNormalization.height);
        var minPoint = new Vector2(pointNormalization.minX, pointNormalization.minY);
        List<Vector2[]> normalizedHoles = new();
        for(int i = 0; i < holes.Count; i++)
        {
            normalizedHoles.Add(NormalizePoints(D, minPoint, holes[i]));
        }

        // Add points
        foreach(var point in pointBinGrid)
            triangleSet.AddPoint(point);

        // Add hole points
        for (int i = 0; i < normalizedHoles.Count; i++)
            foreach(var point in normalizedHoles[i])
                triangleSet.AddPoint(point, i);

        // Fix constrained edges
        // triangleSet.FixConstrainedEdges();

        pointNormalization.normalizedPoints = triangleSet.VertexPositions;
        var vertices = DenormalizePoints(pointNormalization);
        var triangles = triangleSet.Triangles;

        TriangulationResult result;
        result.triangles = triangles.ToArray();
        result.points = vertices;

        return result;
    }

    static PointNormalizationResult NormalizePoints(in Vector2[] pointsToNormalize)
    {
        // Search for min and max values
        float minX, maxX, minY, maxY;

        // Init extreme values
        minX = minY = float.MaxValue;
        maxX = maxY = float.MinValue;

        // Search extreme values
        foreach(var point in pointsToNormalize)
        {
            minX = Mathf.Min(point.x, minX);
            maxX = Mathf.Max(point.x, maxX);
            minY = Mathf.Min(point.y, minY);
            maxY = Mathf.Max(point.y, maxY);
        }

        float width = maxX - minX;
        float height = maxY - minY;

        float D = Mathf.Max(width, height);
        Vector2 minPoint = new Vector2(minX, minY);

        Vector2[] normalizedPoints = NormalizePoints(D, minPoint, pointsToNormalize);

        PointNormalizationResult result = new PointNormalizationResult();
        result.normalizedPoints = normalizedPoints;
        result.minY = minY;
        result.maxY = maxY;
        result.minX = minX;
        result.maxY = maxY;
        result.width = width;
        result.height = height;

        return result;
    }

    static Vector2[] DenormalizePoints(in PointNormalizationResult pointNormalization)
    {
        var result = new Vector2[pointNormalization.normalizedPoints.Length];
        // Denormalize point by point
        var D = Mathf.Max(pointNormalization.width, pointNormalization.height);
        var minPoint = new Vector2(pointNormalization.minX, pointNormalization.minY);

        for (int i = 0; i < result.Length; i++)
            result[i] = minPoint + D * pointNormalization.normalizedPoints[i];

        return result;
    }

    private static Vector2[] NormalizePoints(float D,in Vector2 minPoint, in Vector2[] pointsToNormalize)
    {
        var n = pointsToNormalize.Length;
        var result = new Vector2[n];
        for (int i = 0; i < n; i++)
            result[i] = 1 / D *  (pointsToNormalize[i] - minPoint);

        return result;
    }

    /// <summary>
    /// This function will try to triangulate a polygon specified in COUNTER CLOCKWISE ORDER. This is not garanteed 
    /// to work with polygons specified in clockwise order. This is an implementation of the algorithm in the following
    /// paper to triangulate a PSLG https://cg.cs.uni-bonn.de/backend/v1/files/publications/klein-1996-construction.pdf
    /// </summary>
    /// <param name="polygonIndex">List of indices specifying the polygon</param>
    /// <param name="points">Array of points specifying the actual locatios refered by points in the polygon</param>
    /// <param name="outTriangles">Resulting triangles created by this triangulation</param>
    private static void TriangulatePSLG(LinkedList<int> polygonIndex, in Vector2[] points, ref List<Vector3Int> outTriangles)
    {

        // Check if polygon is already a triangle
        int upTo3Indices = 0;
        foreach(var index in polygonIndex)
        {
            upTo3Indices++;
            if (upTo3Indices > 3)
                break;
        }

        // If three indices, this polygon is already a triangulation
        if (upTo3Indices == 3) 
        {
            int firstIndex = polygonIndex.First.Value;
            int secondIndex = polygonIndex.First.Next.Value;
            int thirdIndex= polygonIndex.First.Next.Next.Value;

            outTriangles.Add(
                new Vector3Int(firstIndex, secondIndex, thirdIndex)
                );
            return;
        }
        if (upTo3Indices < 3) // This is not even a polygon
            return;

        // Now we need to find the edge and the vertex to form the next triangle. We do this by using
        // a function to choose the next edge
        var (start, end) = NextEdge(points, polygonIndex);

        // Now we have to search for every vertex and find one that checks the deluanay condition
        var nextNode = polygonIndex.First;
        while(nextNode != null)
        {
            // Ignore nodes that we already know
            if (nextNode.Value == start.Value || nextNode.Value == end.Value)
            {
                nextNode = nextNode.Next;
                continue;
            }

            Vector2 startToEnd = points[end.Value] - points[start.Value];
            Vector2 startToNext = points[nextNode.Value] - points[start.Value];

            if (!InRightSideOf(startToEnd, startToNext) && 
                IsVisible(polygonIndex, points, nextNode.Value, start.Value) && 
                IsVisible(polygonIndex, points, nextNode.Value, end.Value))
            {
                // Now try to check if this node matches the deluanay condition 
                Vector2 possiblePoint = points[nextNode.Value];
                Vector2 startPoint = points[start.Value];
                Vector2 endPoint = points[end.Value];

                // We have to check if this node also matches the deluanay condition for all other points in the polygon
                bool allDeluanay = true;
                var node = polygonIndex.First;
                while(node != null && allDeluanay)
                {
                    // If this point is part of the original triangle, then just dontinue
                    if (node.Value == start.Value || node.Value == end.Value || node.Value == nextNode.Value)
                    {
                        node = node.Next;
                        continue;
                    }    

                    Vector2 p = points[node.Value];
                    allDeluanay = CheckDeluanayCondition(startPoint, endPoint, possiblePoint, p);
                    node = node.Next;
                }

                // If doesn't match the deluanay condition for every node, return
                if (allDeluanay)
                {
                    // If it does, then we create a new triangle and split our polygon in sub polygons
                    outTriangles.Add(new Vector3Int(start.Value, end.Value, nextNode.Value));
                    var (p0, p1) = SplitPolygonIn(start, end, nextNode);
                    TriangulatePSLG(p0, points, ref outTriangles);
                    TriangulatePSLG(p1, points, ref outTriangles);
                    return;
                }


            }

            nextNode = nextNode.Next;
        }
        Debug.LogError("Couldn't triangulate PSLG properly");

    }

    /// <summary>
    /// Return the next edge to consider represented by its start point index.
    /// </summary>
    /// <param name="points">Array of points specifying a polygon, last edge is points[end] -> points[start] </param>
    /// <param name="start"> Where the specified polygon starts in the array of points </param>
    /// <param name="end"> Where the specified polygon ends in the array of points </param>
    /// <returns> Index of starting end of next edge to choose </returns>
    private static (LinkedListNode<int>, LinkedListNode<int>) NextEdge(in Vector2[] points, LinkedList<int> polygon)
    {
        return (polygon.First, polygon.First.Next);
    }

    /// <summary>
    /// Check if point in `toIndex` is visible from point in `fromIndex`
    /// </summary>
    /// <param name="polygon">Polygon specified as a linked list of vertices</param>
    /// <param name="points">Actual point data</param>
    /// <param name="startIndex">starting point of line of sight</param>
    /// <param name="toIndex">end point of line of sight</param>
    /// <returns></returns>
    private static bool IsVisible(LinkedList<int> polygon, in Vector2[] points, int startIndex,  int toIndex)
    {
        Vector2 from1 = points[startIndex];
        Vector2 to1 = points[toIndex];
        var nextNode = polygon.First;

        // Since we need to check if a point is visible from another, what we really want to know if is there's
        // any line segment in the polygon being intersected by a line coming from and to the target points.
        // This is, there's a line segment interrupting visibility from `start` to `to`
        while(nextNode != null)
        {
            // Ignore if this is one of the elements we already know about
            if (nextNode.Value == startIndex || nextNode.Value == toIndex)
            {
                nextNode = nextNode.Next;
                continue;
            }
            Vector2 from2 = points[nextNode.Value];
            Vector2 to2 = points[(nextNode.Next ?? polygon.First).Value];

            if (LineSegmentsIntersect(from1, to1, from2, to2))
                return false;

            nextNode = nextNode.Next;
        }
        return true;
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
    private static bool CheckDeluanayCondition(in Vector2 v0, in Vector2 v1, in Vector2 v2, in Vector2 p)
    {
        var A = v0;
        var B = v1;
        var C = v2;
        var D = p;

        var mat = new Matrix4x4(
                new Vector4(A.x, B.x, C.x, D.x),
                new Vector4(A.y, B.y, C.y, D.y),
                new Vector4(A.x * A.x + A.y * A.y, B.x * B.x + B.y * B.y, C.x * C.x + C.y * C.y, D.x * D.x + D.y * D.y),
                Vector4.one
            );

        return mat.determinant < 0;
    }

    /// <summary>
    /// Try to split a polygon in two by the specified triangle, the nodes are vertices of this triangles. This generates
    /// Two new polygons
    /// </summary>
    /// <param name="start">Start of triangle, start of an edge in the original polygon</param>
    /// <param name="end">End of triangle, end of an edge in the original polygon</param>
    /// <param name="newVert">new vertex that we choose to form this triangle, usually won't be adjacent to start or end </param>
    /// <returns>Two linked list corresponding to two different polygons formed by this split </returns>
    private static (LinkedList<int>, LinkedList<int>) SplitPolygonIn(LinkedListNode<int> start, LinkedListNode<int> end, LinkedListNode<int> newVert)
    {
        // Create list from end to newVert:
        LinkedList<int> endToNext = new LinkedList<int>();
        var currentNode = end;

        // Creating first list
        while(true)
        {
            endToNext.AddLast(currentNode.Value);
            if (currentNode == newVert)
                break;

            currentNode = currentNode.Next ?? currentNode.List.First;
        }

        // Creating second list, from newVert to start
        LinkedList<int> nextToStart = new LinkedList<int>();
        currentNode = newVert;
        while(true)
        {
            nextToStart.AddLast(currentNode.Value);
            if (currentNode == start)
                break;

            // This line of code help us to traverse this list in a circular fashion
            currentNode = currentNode.Next ?? currentNode.List.First;
        }

        return (endToNext, nextToStart);
    }

    private static bool LineSegmentsIntersect(Vector2 i1, Vector2 f1, Vector2 i2, Vector2 f2)
    {
        var di1 = Cross2D(f1 - i1, i2 - i1);
        var df1 = Cross2D(f1 - i1, f2 - i1);
        var di2 = Cross2D(f2 - i2, i1 - i2);
        var df2 = Cross2D(f2 - i2, f1 - i2);
        var delta = 0.001;

        return ((di1 > 0 == df1 < 0) && (di2 > 0 == df2 < 0))
            || (Mathf.Abs(di1) <= delta && InSegment(i1, f1, i2))
            || (Mathf.Abs(df1) <= delta && InSegment(i1, f1, f2))
            || (Mathf.Abs(di2) <= delta && InSegment(i2, f2, i1))
            || (Mathf.Abs(df2) <= delta && InSegment(i2, f2, f1));
    }

    private static bool InSegment(Vector2 i, Vector2 f, Vector2 p)
    {
        return Mathf.Min(i.x, f.x) < p.x && p.x < Mathf.Max(i.x, f.x) &&
               Mathf.Min(i.y, f.y) < p.y && p.y < Mathf.Max(i.y, f.y);
    }

    private static float Cross2D(Vector2 p1, Vector2 p2)
    {
        return p1.x * p2.y - p2.x * p1.y;
    }

    /// <summary>
    /// Check if `vec` is pointing to the right side of the reference vector `refVec`
    /// </summary>
    /// <param name="refVec">Reference vector from which we want to check if `vec` if in the right side</param>
    /// <param name="vec">Vector we want to check if is in the right side of `refVec`</param>
    /// <returns></returns>
    private static bool InRightSideOf(Vector2 refVec, Vector2 vec) => Cross2D(vec, refVec) > 0;
}
