/*
    Use this classs to split a mesh into multiple meshes by specifying a plane.
    The process can be split in the following steps:
    1. Mesh intersection computation: Mesh intersection points against the mesh in this object
        is computed, outputing the intersection points and triangles.
    2. A new mesh is recomputed: A new mesh is generated following the intersection points to generate
        new borders.
    3. The new mesh is split into as many little meshes is neccesary and new game objects are created with those 
    little meshes.
*/
using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;


[RequireComponent(typeof(MeshFilter))]
public class CutableMeshComponent : MonoBehaviour
{
    public struct TriangleIntersection
    {
        // `piEdge` is the edge being intersected by positioni, an edge is specified by the starting vertex index. If 
        // we have the triangle 0 1 2, and we intersect it in edge 0 -> 1 with position1, then p1Edge == 0
        public Vector3[] position;
        public int triangleIndex;
        public VertexAttributes[] pAttrs;
        public int[] edges; 
    }
    
    public struct VertexAttributes
    {
        public VertexAttributes(Vector2 uvs, Vector3 normal, Vector3 tangent)
        {
            this.uvs = uvs;
            this.normal = normal;
            this.tangent = tangent;
        }

        public Vector2 uvs;
        public Vector3 normal;
        public Vector4 tangent;
    }

    private MeshFilter _meshFilter;

    /// <summary>
    /// Structure used to control adjacency of triangles in mesh. Useful to channel mesh intersection with plane
    /// </summary>
    private int[] _triangleAdjacency;

    // Start is called before the first frame update
    void Start()
    {
        _meshFilter = GetComponent<MeshFilter>();
        var mesh = _meshFilter.mesh;
        _triangleAdjacency = CreateTriangleAdjacencyArray(mesh.triangles);
    }

    // Update is called once per frame
    void Update()
    {
        if (Input.GetKeyDown(KeyCode.Space))
            SplitMesh(new Plane(Vector3.up, 0));
    }

    private void OnDrawGizmos()
    {
        Gizmos.color = Color.red;

        if (_meshFilter) // On editor time _meshFilter is not already set up
        {
            var mesh = _meshFilter.mesh;
            List<TriangleIntersection> intersections;

            // Draw intersection points
            IntersectPlaneToMesh(mesh, new Plane(Vector3.up, 0), out intersections);
            foreach(var intersection in intersections)
            {
                Gizmos.DrawWireSphere(intersection.position[0], 0.01f);
                Gizmos.DrawWireSphere(intersection.position[1], 0.01f);
            }
        }

        // Draw Square simulating plane
        const float sideLen = 10.0f;
        var topLeft = sideLen / 2.0f * Vector3.forward - sideLen / 2.0f * Vector3.right;
        var topRight = sideLen / 2.0f * Vector3.forward + sideLen / 2.0f * Vector3.right;
        var bottomLeft = - sideLen / 2.0f * Vector3.forward - sideLen / 2.0f * Vector3.right;
        var bottomRight = - sideLen / 2.0f * Vector3.forward + sideLen / 2.0f * Vector3.right;

        Gizmos.DrawLine(topLeft, topRight);
        Gizmos.DrawLine(topRight, bottomRight);
        Gizmos.DrawLine(bottomRight, bottomLeft);
        Gizmos.DrawLine(bottomLeft, topLeft);
    }

    /// <summary>
    /// Split this cutable mesh in two meshes by the specified plane*
    /// </summary>
    /// <param name="plane"></param>
    /// <returns></returns>
    private List<Mesh> SplitMesh(Plane plane)
    {
        List<Mesh> result = new List<Mesh>();

        var mesh = _meshFilter.mesh;

        List<TriangleIntersection> intersections;
        if (!IntersectPlaneToMesh(mesh, plane, out intersections)) // No intersection
            return result;

        // Resize array once to avoid multiple array reallocations
        var vertices = mesh.vertices;
        var triangles = mesh.triangles;
        var uvs = mesh.uv;
        var normals = mesh.normals;
        var tangents = mesh.tangents;

        var nextVertexIndex = vertices.Length;
        var nextTrianglesIndex = triangles.Length;

        Array.Resize(ref vertices, vertices.Length + 4 * intersections.Count);
        Array.Resize(ref uvs, uvs.Length + 4 * intersections.Count);
        Array.Resize(ref normals, normals.Length + 4 * intersections.Count); 
        Array.Resize(ref tangents, tangents.Length + 4 * intersections.Count); 

        // Each intersection generates 2 more triangles and modifies one triangle in the previous mesh,
        // so you add 2 * 3 * intersections.Count 
        Array.Resize(ref triangles,  triangles.Length + 3 * 2 * intersections.Count);

        // Use this dictionary to avoid duplication of points that are already created by other triangle intersections.
        // This is a map: Point of intersection -> (side a index, side b index, possible next point)
        var pointGraph = new Dictionary<Vector3, (int, int, HashSet<Vector3>)>();
        int i = 0;
        foreach(var intersection in intersections)
        {
            AddTrianglesFromIntersectingTriangle(intersection, ref vertices, ref uvs, ref normals, ref tangents, ref triangles, nextVertexIndex, nextTrianglesIndex, i++, pointGraph);
        }

        Array.Resize(ref vertices, nextVertexIndex + 2 * pointGraph.Count);
        Array.Resize(ref uvs, nextVertexIndex + 2 * pointGraph.Count);
        Array.Resize(ref normals, nextVertexIndex + 2 * pointGraph.Count);
        Array.Resize(ref tangents, nextVertexIndex + 2 * pointGraph.Count);


        // Now that we added new segments, we can split using disjoints sets,
        // so we know which vertex corresponds to each new mesh
        var (components, nComponents) = ConnectedComponents(triangles, vertices.Length);

        // -- DEBUG ONLY, DELETE LATER -----------------------
        mesh.vertices = vertices;
        mesh.uv = uvs;
        mesh.normals = normals;
        mesh.tangents = tangents;
        mesh.triangles = triangles;
        // ---------------------------------------------------

        // Now that we have all the intersections and recently created points, we will triangulize the set of points.
        // In order to do so, we will create a polygon with this vertices
        // TODO: we assume that there's just a single polygon for now, but there might be more. Think about
        // some ways of cutting an H letter, or slicing a torus

        var polygons = GetPolygonsFromGraph(pointGraph);

        // Now that we have the polygon, we have to convert it to 2D to be triangulated by our algorithm.
        // Note that we don't need to do any back-conversion later, since we only care about indices and not
        // about actual vertices, we already have those.

        // To convert it to 2D, we will create a coordinate frame from plane points, and re-write these
        // points according to that local coordinate frame
        // We know that if we create a coordinate frame i,j,k for this points, we can rewrite each point
        // From the origin, say P' = P - Origin. Then P' = x * i + y * j + z *k, which is the same as
        // [i j k] [x y z]' = P, and finally [x y z]' = [i j k]^-1 x P
        // var origin = polygon[0];
        // for(i = 0; i < polygon.Length; i++)
        // {
        //     var P = polygon[i] - origin;
        // }

        return result;
    }

    /// <summary>
    /// This function will add the input triangle to the specified array of vertices, uvs, normals and tangents
    /// assuming that there's enough space for them. This function will also build the graph formed by the multiple possible
    /// holes generated by the cut, this is necessary so that we can later fill those holes with more triangles.
    /// </summary>
    /// <param name="triangleIntersection"></param>
    /// <param name="vertices"></param>
    /// <param name="uvs"></param>
    /// <param name="normals"></param>
    /// <param name="tangents"></param>
    /// <param name="triangles"></param>
    /// <param name="vertexIndexStart"></param>
    /// <param name="triangleIndexStart"></param>
    /// <param name="intersectionIndex"></param>
    /// <param name="pointGraph"></param>
    /// <param name="pointGraphInverse"></param>
    private void AddTrianglesFromIntersectingTriangle(
        in TriangleIntersection triangleIntersection,
        ref Vector3[] vertices,
        ref Vector2[] uvs,
        ref Vector3[] normals,
        ref Vector4[] tangents,
        ref int[] triangles,
        int vertexIndexStart,
        int triangleIndexStart,
        int intersectionIndex,
        Dictionary<Vector3, (int, int, HashSet<Vector3>)> pointGraph
        )
    {
        // TODO this might be a good spot to register points to project and generate new tessellation to close mesh

        // Add new vertices to vertice array. We need to add four vertices, two per side of splitted mesh, so the resulting meshes
        // Will be disjoint
        var nextVertexIndex = vertexIndexStart + 2 * pointGraph.Count;
        int p1Side1, p1Side2;
        var p1Local = WorldToLocal(triangleIntersection.position[0]);
        var p2Local = WorldToLocal(triangleIntersection.position[1]);

        if (pointGraph.ContainsKey(p1Local)) // TODO: We have to change this to use something with exact precision, this is a POC for now
        {
            HashSet<Vector3> adjacencySet;
            (p1Side1, p1Side2, adjacencySet) = pointGraph[p1Local];
            adjacencySet.Add(p2Local);
        }
        else
        {
            (p1Side1, p1Side2) = (nextVertexIndex, nextVertexIndex + 1);
            var adjacencySet = new HashSet<Vector3>();
            adjacencySet.Add(p2Local);
            pointGraph[p1Local] = (p1Side1, p1Side2, adjacencySet);
            nextVertexIndex += 2;
        }

        int p2Side1, p2Side2;
        if (pointGraph.ContainsKey(p2Local)) // TODO: We have to change this to use something with exact precision, this is a POC for now
        {
            HashSet<Vector3> adjacencySet;
            (p2Side1, p2Side2, adjacencySet) = pointGraph[p2Local];
            adjacencySet.Add(p1Local);
        }
        else
        {
            var adjacencySet = new HashSet<Vector3>();
            adjacencySet.Add(p1Local);
            (p2Side1, p2Side2) = (nextVertexIndex, nextVertexIndex + 1);
            pointGraph[p2Local] = (p2Side1, p2Side2, adjacencySet);
        }



        // This array will map from intersection point (0,1) and side (0, 1) to its corresponding index in the vertex array.
        // For example, A[i,j] is index of point i in side j
        int[,] pointAndSideToIndex = new int[2, 2] { { p1Side1,p1Side2} , { p2Side1, p2Side2} };

        vertices[p1Side1]   = p1Local; // Remember that intersections are computed in world coordinates
        normals[p1Side1]    = triangleIntersection.pAttrs[0].normal;
        uvs[p1Side1]        = triangleIntersection.pAttrs[0].uvs;
        normals[p1Side1]    = triangleIntersection.pAttrs[0].normal;
        tangents[p1Side1]   = triangleIntersection.pAttrs[0].tangent;

        vertices[p2Side1]   = p2Local;
        normals[p2Side1]    = triangleIntersection.pAttrs[1].normal;
        uvs[p2Side1]        = triangleIntersection.pAttrs[1].uvs;
        normals[p2Side1]    = triangleIntersection.pAttrs[1].normal;
        tangents[p2Side1]   = triangleIntersection.pAttrs[1].tangent;

        vertices[p1Side2] = p1Local; // Remember that intersections are computed in world coordinates
        normals[p1Side2] = triangleIntersection.pAttrs[0].normal;
        uvs[p1Side2] = triangleIntersection.pAttrs[0].uvs;
        normals[p1Side2] = triangleIntersection.pAttrs[0].normal;
        tangents[p1Side2] = triangleIntersection.pAttrs[0].tangent;

        vertices[p2Side2] = p2Local;
        normals[p2Side2] = triangleIntersection.pAttrs[1].normal;
        uvs[p2Side2] = triangleIntersection.pAttrs[1].uvs;
        normals[p2Side2] = triangleIntersection.pAttrs[1].normal;
        tangents[p2Side2] = triangleIntersection.pAttrs[1].tangent;

        var v0Index = triangles[triangleIntersection.triangleIndex];
        var v1Index = triangles[triangleIntersection.triangleIndex+1];
        var v2Index = triangles[triangleIntersection.triangleIndex+2];

        // Check which vertex of previous triangle was alone in the other side of the plane
        int aloneIndex = -1; // This is the vertex isolated in one side of the plane
        int midPrevIndex = -1;  // This is the index of the intersection position that's in the edge of aloneIndex - 1 -> aloneIndex
        int midNextIndex = -1;  // This is the index of the intersection position that's in the edge of aloneIndex -> aloneIndex + 1
        for (int vert = 0; vert < 3; vert ++)
        {
            if ((triangleIntersection.edges[0] + 1) % 3 == vert && triangleIntersection.edges[1] == vert)
            {
                aloneIndex = vert;
                midPrevIndex = 0;
                midNextIndex = 1;
                break;
            }
            else if ((triangleIntersection.edges[1] + 1) % 3 == vert && triangleIntersection.edges[0] == vert)
            {
                aloneIndex = vert;
                midPrevIndex = 1;
                midNextIndex = 0;
                break;
            }
        }
        int prevIndex = aloneIndex - 1 < 0 ? 2 : aloneIndex - 1;
        int nextIndex = (aloneIndex + 1) % 3;
        Debug.Assert(aloneIndex != -1, "Couldn't find index alone in one side of plane");

        int[] originalTriangle = new int[] { triangles[triangleIntersection.triangleIndex], triangles[triangleIntersection.triangleIndex + 1], triangles[triangleIntersection.triangleIndex + 2] };

        // Update old triangle
        // Leave first as is:
        // triangles[triangleIntersection.triangleIndex + aloneIndex] = triangles[triangleIntersection.triangleIndex + aloneIndex]; 
        triangles[triangleIntersection.triangleIndex + prevIndex] = pointAndSideToIndex[midPrevIndex, 0];
        triangles[triangleIntersection.triangleIndex + nextIndex] = pointAndSideToIndex[midNextIndex, 0];

        // Create two new triangles
        //      One triangle
        var nextTriangleIndex = triangleIndexStart + 6 * intersectionIndex;
        triangles[nextTriangleIndex + aloneIndex] = pointAndSideToIndex[midNextIndex, 1];
        triangles[nextTriangleIndex + prevIndex] = pointAndSideToIndex[midPrevIndex, 1];
        triangles[nextTriangleIndex + nextIndex] = originalTriangle[nextIndex];
        nextTriangleIndex += 3;

        //      Another triangle
        triangles[nextTriangleIndex + prevIndex] = originalTriangle[prevIndex];
        triangles[nextTriangleIndex + nextIndex] = originalTriangle[nextIndex];
        triangles[nextTriangleIndex + aloneIndex] = pointAndSideToIndex[midPrevIndex, 1];
    }

    bool IntersectPlaneToMesh(in Mesh mesh, in Plane plane, out List<TriangleIntersection> intersectingTriangles)
    {
        intersectingTriangles = new List<TriangleIntersection>();

        var meshBounds = mesh.bounds;

        // Bounding box points to check if they're all in the same side of plane
        var boundPoint1 = LocalToWorld(meshBounds.min);
        var boundPoint2 = LocalToWorld(meshBounds.max);
        var boundPoint3 = new Vector3(boundPoint1.x, boundPoint1.y, boundPoint2.z);
        var boundPoint4 = new Vector3(boundPoint1.x, boundPoint2.y, boundPoint1.z);
        var boundPoint5 = new Vector3(boundPoint2.x, boundPoint1.y, boundPoint1.z);
        var boundPoint6 = new Vector3(boundPoint1.x, boundPoint2.y, boundPoint2.z);
        var boundPoint7 = new Vector3(boundPoint2.x, boundPoint1.y, boundPoint2.z);
        var boundPoint8 = new Vector3(boundPoint2.x, boundPoint2.y, boundPoint1.z);

        Vector3[] boundPoints = {
            boundPoint1, boundPoint2, boundPoint3, boundPoint4,
            boundPoint5, boundPoint6, boundPoint7, boundPoint8
        };

        // If both extents of mesh are in the same side of the plane, the plane does not intersects
        var firstBoundPtsSign = SideOfPlane(plane.normal, plane.distance, boundPoint1);

        bool allSameSide = true;
        for (uint i = 1; i < boundPoints.Length && allSameSide; i++)
        {
            allSameSide = firstBoundPtsSign == SideOfPlane(plane.normal, plane.distance, boundPoints[i]);
        }

        if (allSameSide)
            return false;

        // Perform actual intersection now that we know that we actually might be intersecting this mesh
        var triangles = mesh.triangles;
        var normals = mesh.normals;
        var tangents = mesh.tangents;
        var uvs = mesh.uv;
        var vertices = mesh.vertices;
        for (int i = 0; i < mesh.triangles.Length; i += 3 )
        {
            TriangleIntersection result;
            if (IntersectPlaneToTrianglev2(triangles, vertices, normals, tangents, uvs, i, plane, out result))
            {
                result.triangleIndex = i;
                intersectingTriangles.Add(result);
            }
        }

        return intersectingTriangles.Count != 0;
    }

    private bool IntersectPlaneToTrianglev2(
        in int[] triangles,
        in Vector3[] vertices,
        in Vector3[] normals,
        in Vector4[] tangents,
        in Vector2[] uvs,
        int triangleIndex,
        in Plane plane,
        out TriangleIntersection intersection
        )
    {
        var worldVertices = new Vector3[] { LocalToWorld(vertices[triangles[triangleIndex]]), LocalToWorld(vertices[triangles[triangleIndex + 1]]), LocalToWorld(vertices[triangles[triangleIndex + 2]]) };

        // if all points lay in the same side or some of them touch the plane, this triangle is not intersecting the plane
        var v0PlaneSide = SideOfPlane(plane.normal, plane.distance, worldVertices[0]);
        var v1PlaneSide = SideOfPlane(plane.normal, plane.distance, worldVertices[1]);
        var v2PlaneSide = SideOfPlane(plane.normal, plane.distance, worldVertices[2]);
        intersection = new TriangleIntersection();
        if ((v0PlaneSide == v1PlaneSide && v1PlaneSide == v2PlaneSide) || v0PlaneSide == 0 || v1PlaneSide == 0 || v2PlaneSide == 0)
            return false;

        // Initialize result
        intersection = new TriangleIntersection();
        intersection.position = new Vector3[2];
        intersection.pAttrs = new VertexAttributes[2];
        intersection.triangleIndex = triangleIndex;
        intersection.edges = new int[2];

        // Iterate over edges of this triangle to try to intersect this plane
        var nextIntersection = 0;
        for (var edge = 0; edge < 3 && nextIntersection < 2; edge ++)
        {
            var vStartIndex = triangles[triangleIndex + edge];
            var vEndIndex = triangles[triangleIndex + (edge + 1) % 3];
            var vStart = worldVertices[edge];
            var vEnd = worldVertices[(edge+1) % 3];
            Vector3 intersectionPoint;
            float t;
            if (IntersectSegmentToPlane(vStart, vEnd, plane.normal, plane.distance, out intersectionPoint, out t) && 0 < t && t < 1)
            {
                intersection.position[nextIntersection] = intersectionPoint;
                intersection.edges[nextIntersection] = edge;
                intersection.pAttrs[nextIntersection] = InterpolateVertexAttributes(new VertexAttributes(uvs[vStartIndex], normals[vStartIndex], tangents[vStartIndex]), 
                                                                                    new VertexAttributes(uvs[vEndIndex], normals[vEndIndex], tangents[vEndIndex]), t);
                nextIntersection++;
            }
        }

        Debug.Assert(nextIntersection == 2, "Invalid amount of intersections in a triangle, should be 2");

        return true;
    }

    private bool IntersectSegmentToPlane(in Vector3 v0, in Vector3 v1, in Vector3 planeNormal, float planeDistance, out Vector3 intersectionPoint, out float t)
    {
        float d0 = DistancePointToPlane(v0, planeNormal, planeDistance);
        float d1 = DistancePointToPlane(v1, planeNormal, planeDistance);

        if (d0 * d1 >= 0 ) // Points in the same side of plane
        {
            intersectionPoint = Vector3.zero;
            t = 0;
            return false;
        }

        t = d0 / (d0 - d1);
        intersectionPoint = v0 + t * (v1 - v0);
        
        return true;
    }

    private float DistancePointToPlane(in Vector3 v0, in Vector3 planeNormal, in float planeDistance)
    {
        return Vector3.Dot(v0, planeNormal) + planeDistance;
    }

    private float SideOfPlane(in Vector3 planeNormal, float planeDist,  in Vector3 point)
    {

        var planeEval = Vector4.Dot(
                new Vector4(point.x, point.y, point.z, planeDist),
                new Vector4(planeNormal.x, planeNormal.y, planeNormal.z, 1)
            );
        return Mathf.Abs(planeEval) < 0.0001 ? 0 : Mathf.Sign(planeEval);
    }

    /// <summary>
    /// Transform a point from local model coordinates to world coordinates
    /// </summary>
    /// <param name="point"> Point in model coordinates to transform to world coordinates </param>
    /// <returns> Point transformed from local model coordinates to world coordinates </returns>
    private Vector3 LocalToWorld(in Vector3 point)
    {
        // Create point in homogeneous coordinates
        Vector4 pointHomoCoords;
        pointHomoCoords.x = point.x;
        pointHomoCoords.y = point.y;
        pointHomoCoords.z = point.z;
        pointHomoCoords.w = 1.0f;

        // Transform point to world coordinates
        pointHomoCoords = transform.localToWorldMatrix * pointHomoCoords;

        // Transform back to vector3
        Vector3 result;
        result.x = pointHomoCoords.x;
        result.y = pointHomoCoords.y;
        result.z = pointHomoCoords.z;

        return result;
    }

    /// <summary>
    /// Transform a point from woolrd coordinates to local model coordinates
    /// </summary>
    /// <param name="point"> Point in world coordinates to transform to model coordinates </param>
    /// <returns> Point transformed from world coordinates to model coordinates </returns>
    private Vector3 WorldToLocal(in Vector3 point)
    {
        // Create point in homogeneous coordinates
        Vector4 pointHomoCoords;
        pointHomoCoords.x = point.x;
        pointHomoCoords.y = point.y;
        pointHomoCoords.z = point.z;
        pointHomoCoords.w = 1.0f;

        // Transform point to world coordinates
        pointHomoCoords = transform.worldToLocalMatrix * pointHomoCoords;

        // Transform back to vector3
        Vector3 result;
        result.x = pointHomoCoords.x;
        result.y = pointHomoCoords.y;
        result.z = pointHomoCoords.z;

        return result;
    }

    /// <summary>
    /// Compute triangle normal according to unity's winding order and coordinates system   
    /// </summary>
    /// <param name="v0">first triangle vertex</param>
    /// <param name="v1">second triangle vertex</param>
    /// <param name="v2">third triangle vertex</param>
    /// <param name="clockwise">If winding order is clockwise or not. Unity's winding order is clockwise </param>
    /// <returns>Normal of  the specified triangle </returns>
    public Vector3 TriangleNormal(Vector3 v0, Vector3 v1, Vector3 v2, bool clockwise = true)
    {
        Vector3 u = v1 - v0;    // edge v0 -> v1
        Vector3 v = v2 - v0;    // edge v0 -> v2

        // Unity uses Clockwise winding order to determine front-facing triangles
        // Unity uses a left-handed coordinate system
        // the normal faces front
        Vector3 normal = (clockwise ? 1 : -1) * Vector3.Cross(u, v).normalized;
        return normal;
    }

    public VertexAttributes InterpolateVertexAttributes(VertexAttributes v0Attr, VertexAttributes v1Attr, float t)
    {
        if (!(0 <= t && t <= 1))
            Debug.Assert(false, "Bad intersection");

        Debug.Assert(0 <= t && t <= 1, $"Inconsistent interpolation value: {t}");
        VertexAttributes newAttrs;
        newAttrs.normal = Vector3.Lerp(v0Attr.normal, v1Attr.normal, t);
        newAttrs.tangent = Vector4.Lerp(v0Attr.tangent, v1Attr.tangent, t);
        newAttrs.uvs = Vector2.Lerp(v0Attr.uvs, v1Attr.uvs, t);

        return newAttrs;
    }


    private struct DisjointSet
    {
        public int[] parents;
        public int nVertices;

        public DisjointSet(int nVertices)
        {
            parents = new int[nVertices];
            this.nVertices = nVertices;

            // Init array of parents as your own parent
            for (int i = 0; i < nVertices; i++)
            {
                parents[i] = i;
            }
        }

        /// <summary>
        /// Get root of node i, which must be in range [0, nVertices)
        /// </summary>
        /// <param name="i">node whose parent you want</param>
        public int Root(int i)
        {
            if (parents[i] == i)
                return i;

            // If you're not your own parent, search your parent's parent
            var actualParent = Root(parents[i]);

            // Optimization to reduce height of resulting parent tree 
            parents[i] = actualParent;

            return actualParent;
        }

        /// <summary>
        /// Connect to subsets represented by their members i and j
        /// </summary>
        /// <param name="i">Member of subset 1</param>
        /// <param name="j">Member of subset 2</param>
        public void Merge(int i, int j)
        {
            // Make i be parent of j's parent
            var jRoot = Root(j);
            var iRoot = Root(i);
            if (jRoot == iRoot)
                return; // nothing to do, already in the same set

            parents[jRoot] = iRoot;
        }
    }

    /// <summary>
    /// A disjoint sets based connected components algorithm to separate a triangle-based graph into multiple components
    /// </summary>
    /// <param name="triangles">array of triangles, each 3 ints is a triangle</param>
    /// <param name="nVertices">ammount of vertices, this is the size of the returning array</param>
    /// <returns>Array of size `nVertices`, each index corresponds to the connected component for this vertex, and the amount of components found</returns>
    private (int[], int) ConnectedComponents(in int [] triangles, int nVertices)
    {
        DisjointSet disjointSet = new DisjointSet(nVertices);
        for (int i = 0; i < triangles.Length; i += 3)
        {
            var v1 = triangles[i];
            var v2 = triangles[i+1];
            var v3 = triangles[i+2];

            disjointSet.Merge(v1, v2);
            disjointSet.Merge(v1, v3);
        }

        SortedSet<int> components = new SortedSet<int>();

        for(int i = 0; i < nVertices; i++)
        {
            var iRoot = disjointSet.Root(i); // flatten entire array of parents
            components.Add(iRoot);
        }

        return (disjointSet.parents, components.Count);
    }

    /// <summary>
    /// Create a triangle adjacency array from a verex array. Each position i will have 
    /// the index of triangle adjacent to edge vertices[i] -> vertices[(i+1) % 3]
    /// </summary>
    /// <param name="vertices">List of triangles specified as soup of triangles</param>
    /// <returns>Adjacency array for each triangle</returns>
    private static int[] CreateTriangleAdjacencyArray(int[] vertices)
    {
        int[] adjacency = new int[vertices.Length];

        // Fill everything with -1 as default value
        Array.Fill<int>(adjacency, -1);

        // For every triangle in array...
        for (int i = 0; i < adjacency.Length; i+=3)
        {
            // If we alreadoy found 3 triangles, then we have nothing more to look for
            int count = 0;

            // For every following triangle while we have neighbors to find...
            for (int j = i+3; j < adjacency.Length && count < 3; j+=3)
            {
                // For each edge in triangle i...
                for (int i_edge = 0; i_edge < 3; i_edge++)
                {
                    var i_start = vertices[i + i_edge];
                    var i_end = vertices[i + (i_edge + 1) % 3];

                    // For each edge in triangle j...
                    for (int j_edge = 0; j_edge < 3; j_edge++)
                    {
                        var j_start = vertices[j];
                        var j_end = vertices[j + (j_edge + 1) % 3];

                        // If we found a matching edge, save it and look for the next edge
                        if ((j_start == i_start && j_end == i_end) || (j_end == i_start && j_start == i_end))
                        {
                            adjacency[i + i_edge] = j;
                            count++;
                            break;
                        }
                    }
                }
                
            }
        }

        return adjacency;
    }

    /// <summary>
    /// Generate the linked lists corresponding to polygons in the graph
    /// </summary>
    /// <param name="pointGraph"></param>
    /// <returns></returns>
    private static List<LinkedList<(int,int)>> GetPolygonsFromGraph(Dictionary<Vector3, (int, int, HashSet<Vector3>)> pointGraph)
    {
        HashSet<Vector3> memo = new HashSet<Vector3>();
        List<LinkedList<(int, int)>> result = new List<LinkedList<(int, int)>>();

        foreach(var v in pointGraph.Keys)
        {
            // Check if already visited this polygon
            if (memo.Contains(v))
                continue;

            // Init next polygon to create
            var nextPoly = new LinkedList<(int, int)>();
            Vector3 nextNode = v;
            int sideA, sideB;
            HashSet<Vector3> nextNodeNbrs;
            (sideA, sideB, nextNodeNbrs) = pointGraph[nextNode];

            while(!memo.Contains(nextNode))
            {
                // Add this node to polygon
                memo.Add(nextNode);
                nextPoly.AddLast((sideA, sideB));

                // Search for next node
                foreach(var nbr in nextNodeNbrs)
                {
                    if (memo.Contains(nbr))
                        continue;

                    nextNode = nbr;
                    (sideA, sideB, nextNodeNbrs) = pointGraph[nbr];
                    break;
                }
            }


            // TODO: Should check if this is a closed polygon, 
            // you can do it by just checking if first and last of nextPoly linked list contain each other
            if (nextPoly.Count > 2) // if an actual polygon, might have less vertices in cases where the shape is not closed 
                result.Add(nextPoly);
        }

        return result;
    }
}
