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

    public static TriangulationResult Triangulate(in Vector2[] pointsToTriangulate, in List<Vector2[]> holes)
    {
        var pointNormalization = NormalizePoints(pointsToTriangulate);
        var pointBinGrid = new PointBinGrid(pointNormalization.normalizedPoints);
        var triangleSet = new DeluanayTriangleSet((uint) pointsToTriangulate.Length, holes.Count); 

        // Normalize holes using the same normalization process for the input points
        var D = Mathf.Max(pointNormalization.width, pointNormalization.height);
        var minPoint = new Vector2(pointNormalization.minX, pointNormalization.minY);
        List<Vector2[]> normalizedHoles = new(holes.Count);
        for(int i = 0; i < normalizedHoles.Count; i++)
        {
            normalizedHoles[i] = NormalizePoints(D, minPoint, holes[i]);
        }

        foreach(var point in pointBinGrid)
            triangleSet.AddPoint(point);

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
}
