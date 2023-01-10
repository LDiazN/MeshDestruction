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
        public Vector3 position1; // always intersect triangle in at the least two positions
        public Vector3 position2;
        public uint v0Index;         // Index of vertex 0 in vertices array
        public uint v1Index;         // Index of vertex 1 in vertices array
        public uint v2Index;         // Index of vertex 2 in vertices array
        public uint v0TriangleIndex; // index in triangles array correspoing to the index of v0
        public uint v1TriangleIndex; // index in triangles array correspoing to the index of v1
        public uint v2TriangleIndex; // index in triangles array correspoing to the index of v
        public VertexAttributes p1Attrs;
        public VertexAttributes p2Attrs;
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
            IntersectPlaneToMesh(mesh, Vector3.up, 0, out intersections);
            foreach(var intersection in intersections)
            {
                Gizmos.DrawWireSphere(intersection.position1, 0.01f);
                Gizmos.DrawWireSphere(intersection.position2, 0.01f);
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
        if (!IntersectPlaneToMesh(mesh, plane.normal, plane.distance, out intersections)) // No intersection
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
        Dictionary<Vector3, (int, int, Vector3?)> existingPoints = new Dictionary<Vector3, (int, int, Vector3?)>();
        int i = 0;
        foreach(var intersection in intersections)
        {
            AddTrianglesFromIntersectingTriangle(intersection, plane, ref vertices, ref uvs, ref normals, ref tangents, ref triangles, nextVertexIndex, nextTrianglesIndex, i++, existingPoints);
        }

        Array.Resize(ref vertices, nextVertexIndex + 2 * existingPoints.Count);
        Array.Resize(ref uvs, nextVertexIndex + 2 * existingPoints.Count);
        Array.Resize(ref normals, nextVertexIndex + 2 * existingPoints.Count);
        Array.Resize(ref tangents, nextVertexIndex + 2 * existingPoints.Count);


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
        Vector3[] polygon = new Vector3[existingPoints.Count];
        Vector3? nextVertex = vertices[nextVertexIndex];
        bool isPolygon = nextVertex != null;
        polygon[0] = vertices[nextVertexIndex];
        for(i = 1; i < existingPoints.Count && isPolygon; i++)
        {
            (_, _, nextVertex) = existingPoints[nextVertex ?? Vector3.zero];
            isPolygon = nextVertex != null;
            polygon[i] = nextVertex ?? Vector3.zero;
        }

        if (!isPolygon)
            return result;
        // Now that we have the polygon, we have to convert it to 2D to be triangulated by our algorithm.
        // Note that we don't need to do any back-conversion later, since we only care about indices and not
        // about actual vertices, we already have those.

        // To convert it to 2D, we will create a coordinate frame from plane points, and re-write these
        // points according to that local coordinate frame
        // We know that if we create a coordinate frame i,j,k for this points, we can rewrite each point
        // From the origin, say P' = P - Origin. Then P' = x * i + y * j + z *k, which is the same as
        // [i j k] [x y z]' = P, and finally [x y z]' = [i j k]^-1 x P
        var origin = polygon[0];
        for(i = 0; i < polygon.Length; i++)
        {
            var P = polygon[i] - origin;

        }

        return result;
    }

    private void AddTrianglesFromIntersectingTriangle(
        in TriangleIntersection triangleIntersection, 
        in Plane plane,
        ref Vector3[] vertices, 
        ref Vector2[] uvs,
        ref Vector3[] normals,
        ref Vector4[] tangents,
        ref int[] triangles, 
        int vertexIndexStart, 
        int triangleIndexStart, 
        int intersectionIndex,
        Dictionary<Vector3, (int,int, Vector3?)> existingPoints
        )
    {
        // TODO this might be a good spot to register points to project and generate new tessellation to close mesh

        // Add new vertices to vertice array. We need to add four vertices, two per side of splitted mesh, so the resulting meshes
        // Will be disjoint
        var nextVertexIndex = vertexIndexStart + 2 * existingPoints.Count;
        int p1Side1, p1Side2;
        var p1Local = WorldToLocal(triangleIntersection.position1);
        var p2Local = WorldToLocal(triangleIntersection.position2);

        if (existingPoints.ContainsKey(p1Local)) // TODO: We have to change this to use something with exact precision, this is a POC for now
        {
            Vector3? next;
            (p1Side1, p1Side2, next) = existingPoints[p1Local];
            if (next == null)
                existingPoints[p1Local] = (p1Side1, p1Side2, p2Local);
        }
        else
        {
            (p1Side1, p1Side2) = (nextVertexIndex, nextVertexIndex + 1);
            existingPoints[p1Local] = (p1Side1, p1Side2, p2Local);
            nextVertexIndex += 2;
        }

        int p2Side1, p2Side2;
        if (existingPoints.ContainsKey(p2Local)) // TODO: We have to change this to use something with exact precision, this is a POC for now
        {
            (p2Side1, p2Side2, _) = existingPoints[p2Local];
        }
        else
        {
            (p2Side1, p2Side2) = (nextVertexIndex, nextVertexIndex + 1);
            existingPoints[p2Local] = (p2Side1, p2Side2, null);
        }

        vertices[p1Side1]   = p1Local; // Remember that intersections are computed in world coordinates
        normals[p1Side1]    = triangleIntersection.p1Attrs.normal;
        uvs[p1Side1]        = triangleIntersection.p1Attrs.uvs;
        normals[p1Side1]    = triangleIntersection.p1Attrs.normal;
        tangents[p1Side1]   = triangleIntersection.p1Attrs.tangent;

        vertices[p2Side1]   = p2Local;
        normals[p2Side1]    = triangleIntersection.p2Attrs.normal;
        uvs[p2Side1]        = triangleIntersection.p2Attrs.uvs;
        normals[p2Side1]    = triangleIntersection.p2Attrs.normal;
        tangents[p2Side1]   = triangleIntersection.p2Attrs.tangent;

        vertices[p1Side2] = p1Local; // Remember that intersections are computed in world coordinates
        normals[p1Side2] = triangleIntersection.p1Attrs.normal;
        uvs[p1Side2] = triangleIntersection.p1Attrs.uvs;
        normals[p1Side2] = triangleIntersection.p1Attrs.normal;
        tangents[p1Side2] = triangleIntersection.p1Attrs.tangent;

        vertices[p2Side2] = p2Local;
        normals[p2Side2] = triangleIntersection.p2Attrs.normal;
        uvs[p2Side2] = triangleIntersection.p2Attrs.uvs;
        normals[p2Side2] = triangleIntersection.p2Attrs.normal;
        tangents[p2Side2] = triangleIntersection.p2Attrs.tangent;

        // Check which vertex of previous triangle was alone in the other side of the plane
        var v0Side = SideOfPlane(plane.normal, plane.distance, LocalToWorld(vertices[triangleIntersection.v0Index]));
        var v1Side = SideOfPlane(plane.normal, plane.distance, LocalToWorld(vertices[triangleIntersection.v1Index]));
        var v2Side = SideOfPlane(plane.normal, plane.distance, LocalToWorld(vertices[triangleIntersection.v2Index]));

        int aloneIndex = -1, same1Index = -1, same2Index = -1;

        if (v0Side != v1Side && v0Side != v2Side)
        {
            aloneIndex = (int) triangleIntersection.v0Index;

            same1Index = (int)triangleIntersection.v1Index;

            same2Index = (int)triangleIntersection.v2Index;
        }
        else if (v1Side != v0Side && v1Side != v2Side)
        {
            aloneIndex = (int)triangleIntersection.v1Index;

            same1Index = (int)triangleIntersection.v0Index;

            same2Index = (int)triangleIntersection.v2Index;
        }
        else if (v2Side != v0Side && v2Side != v1Side)
        {
            aloneIndex = (int)triangleIntersection.v2Index;

            same1Index = (int)triangleIntersection.v0Index;

            same2Index = (int)triangleIntersection.v1Index;
        }
        else
        {
            Debug.Assert(false, "wtf"); // TODO debug weird bug where sometimes intersected triangles are all in the same side of the plane
        }

        // Compute current triangle normal, so we know where our triangle should face
        var originalNormal = TriangleNormal(vertices[triangleIntersection.v0Index], vertices[triangleIntersection.v1Index], vertices[triangleIntersection.v2Index]);


        // Update old triangle
        var (i, j, k) = ComputeWindingOrder(p2Side1, p1Side1, aloneIndex, originalNormal, vertices);
        triangles[triangleIntersection.v0TriangleIndex] = i;
        triangles[triangleIntersection.v2TriangleIndex] = j;
        triangles[triangleIntersection.v1TriangleIndex] = k;

        // Create two new triangles
        //      One triangle
        var nextTriangleIndex = triangleIndexStart + 6 * intersectionIndex;
        (i, j, k) = ComputeWindingOrder(same2Index, p1Side2, p2Side2, originalNormal, vertices);
        triangles[nextTriangleIndex++] = i;
        triangles[nextTriangleIndex++] = j;
        triangles[nextTriangleIndex++] = k;

        //      Another triangle
        (i, j, k) = ComputeWindingOrder(same2Index, same1Index, p1Side2, originalNormal, vertices);
        triangles[nextTriangleIndex++] = i;
        triangles[nextTriangleIndex++] = j;
        triangles[nextTriangleIndex++] = k;
    }

    bool IntersectPlaneToMesh(in Mesh mesh, in Vector3 planeNormal, float planeDistance, out List<TriangleIntersection> intersectingTriangles)
    {
        Debug.Assert(planeNormal.magnitude == 1, "Plane Normal should be normalized to ensure consistent results");

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
        var firstBoundPtsSign = SideOfPlane(planeNormal, planeDistance, boundPoint1);

        bool allSameSide = true;
        for (uint i = 1; i < boundPoints.Length && allSameSide; i++)
        {
            allSameSide = firstBoundPtsSign == SideOfPlane(planeNormal, planeDistance, boundPoints[i]);
        }

        if (allSameSide)
            return false;

        // Perform actual intersection now that we know that we actually might be intersecting this mesh
        var triangles = mesh.triangles;
        for(uint i = 0; i < mesh.triangles.Length; i += 3 )
        {
            
            var v0 = LocalToWorld(mesh.vertices[triangles[i]]);
            var v2 = LocalToWorld(mesh.vertices[triangles[i + 1]]);
            var v1 = LocalToWorld(mesh.vertices[triangles[i + 2]]);
            VertexAttributes v0Attrs = new VertexAttributes(mesh.uv[triangles[i]], mesh.normals[triangles[i]], mesh.tangents[triangles[i]]);
            VertexAttributes v1Attrs = new VertexAttributes(mesh.uv[triangles[i + 2]], mesh.normals[triangles[i + 2]], mesh.tangents[triangles[i + 2]]);
            VertexAttributes v2Attrs = new VertexAttributes(mesh.uv[triangles[i + 1]], mesh.normals[triangles[i + 1]], mesh.tangents[triangles[i + 1]]);

            TriangleIntersection result;
            if (IntersectPlaneToTriangle(v0, v1, v2, planeNormal, planeDistance, v0Attrs, v1Attrs, v2Attrs, out result))
            {
                result.v0TriangleIndex = i;
                result.v2TriangleIndex = i+1;
                result.v1TriangleIndex = i+2;
                result.v0Index = (uint) triangles[i];
                result.v2Index = (uint) triangles[i+1];
                result.v1Index = (uint) triangles[i+2];
                intersectingTriangles.Add(result);
            }
        }

        return intersectingTriangles.Count != 0;
    }

    private bool IntersectPlaneToTriangle(
        in Vector3 v0, in Vector3 v1, in Vector3 v2, 
        in Vector3 planeNormal, float planeDistance, 
        in VertexAttributes v0Attrs, in VertexAttributes v1Attrs, in VertexAttributes v2Attrs, 
        out TriangleIntersection intersection)
    {
        Vector3 intersectionPoint;
        List<(Vector3, VertexAttributes)> intersections = new List<(Vector3, VertexAttributes)>(); // Should not contain more than 2 points
        float t;

        if (IntersectSegmentToPlane(v0, v1, planeNormal, planeDistance, out intersectionPoint, out t))
        {
            // Interpolate properties of this vertex
            VertexAttributes newAttrs = InterpolateVertexAttributes(v0Attrs, v1Attrs, t);
            intersections.Add((intersectionPoint, newAttrs));

        }

        if (IntersectSegmentToPlane(v0, v2, planeNormal, planeDistance, out intersectionPoint, out t))
        {
            VertexAttributes newAttrs = InterpolateVertexAttributes(v0Attrs, v2Attrs, t);
            intersections.Add((intersectionPoint, newAttrs));
        }

        // Don't check if you already have 2 intersections, save some operations
        if (intersections.Count < 2 && IntersectSegmentToPlane(v1, v2, planeNormal, planeDistance, out intersectionPoint, out t))
        {
            VertexAttributes newAttrs = InterpolateVertexAttributes(v1Attrs, v2Attrs, t);
            intersections.Add((intersectionPoint, newAttrs));
        }

        // Sanity check
        Debug.Assert(intersections.Count == 0 || intersections.Count == 2, "Invalid amount of intersections between ray and triangle");

        intersection = new TriangleIntersection();
        if (intersections.Count == 0)
            return false;

        (intersection.position1, intersection.p1Attrs) = intersections[0];
        (intersection.position2, intersection.p2Attrs) = intersections[1];

        return intersections.Count > 0;
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

        return Mathf.Sign(
            Vector4.Dot(
                new Vector4(point.x, point.y, point.z, planeDist),
                new Vector4(planeNormal.x, planeNormal.y, planeNormal.z, 1)
            )
        );
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
    /// Compute correct winding order such that the resulting triangle has a normal vector that points in the same direction of the provided  
    /// expected normal
    /// </summary>
    /// <param name="v0Index">index of </param>
    /// <param name="v1Index"></param>
    /// <param name="v2Index"></param>
    /// <param name="expectedNormal"></param>
    /// <param name="vertices"></param>
    /// <returns></returns>
    private (int,int,int) ComputeWindingOrder(int v0Index, int v1Index, int v2Index, in Vector3 expectedNormal, in Vector3[] vertices)
    {
        var v0 = vertices[v0Index];
        var v1 = vertices[v1Index];
        var v2 = vertices[v2Index];

        var triangleNormal = TriangleNormal(v0, v1, v2);
        if (Vector3.Dot(triangleNormal, expectedNormal) < 0) // Pointing in same direction
            return (v0Index, v1Index, v2Index);
        else // Flip last two in order to change order
            return (v0Index, v2Index, v1Index);
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
            for (int j = i+1; j < adjacency.Length && count < 3; j+=3)
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
}
