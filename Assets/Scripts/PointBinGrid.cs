using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;

/// <summary>
/// Implements a data structure to store and retrieve points easier
/// and efficient, used to compute deluanay triangulation and constrained triangulation
/// </summary>
public class PointBinGrid : IEnumerable<Vector3>
{
    // Number of points to store
    private uint _nPoints;

    // Return M from MxM, the size of bin grid
    private uint _gridSize;

    // Bin grid
    private List<Vector3>[] _bins;

    // bottom left corner of given points
    private Vector2 _bottomLeft;

    // top right corner of given points
    private Vector2 _topRight;

    // Get number of points
    public uint NPoints { get { return _nPoints; } }
    public uint GridSize { get { return _gridSize; } } 

    // Get Width of grid
    public float GridWidth { get { return _topRight.x - _bottomLeft.x; } }
    public float GridHeight { get { return _topRight.x - _bottomLeft.x; } }

    // List of currently stored bins
    public List<Vector3>[] Bins { get { return _bins; } }

    public PointBinGrid(Vector2[] points)
    {
        _nPoints = (uint) points.Length;
        _gridSize = (uint) Mathf.Pow(NPoints, 1.0f/4.0f);
        _bins = new List<Vector3>[_gridSize * _gridSize];

        // Init internal lists
        for(uint i = 0; i < _bins.Length; i++)
            _bins[i] = new List<Vector3>();

        (_bottomLeft, _topRight) = CalculateBounds(points);

        foreach (var point in points)
        {
            var binIndex = BinIndex(point);
            _bins[binIndex].Add(point);
        }
    }

    /// <summary>
    /// Return extents for the provided array of points
    /// </summary>
    /// <param name="points"> Array of points whose extents you want to compute </param>
    /// <returns>(bottom left, top right) extents of a rectangle covering all points</returns>
    static private (Vector2, Vector2) CalculateBounds(Vector2[] points)
    {
        Debug.Assert(points.Length > 1, "There should be at the least two points");
        Vector3 bottomLeft = points[0], topRight = points[0];

        foreach(var point in points)
        {
            bottomLeft.x = Mathf.Min(point.x, bottomLeft.x);
            bottomLeft.y = Mathf.Min(point.y, bottomLeft.y);

            topRight.x = Mathf.Max(point.x, topRight.x);
            topRight.y = Mathf.Max(point.y, topRight.y);
        }

        return (bottomLeft, topRight);
    }

    /// <summary>
    /// Compute row and column of provided point in the grid
    /// </summary>
    /// <param name="position">Position in world coordinates</param>
    /// <returns> (row, column) in the grid for the specified position</returns>
    private (uint, uint) ComputeGridPosition(Vector2 position)
    {
        // Convert current position to grid space
        var posGridSpace = position - _bottomLeft;
        var row = 0.99f * _gridSize * posGridSpace.y / GridHeight;
        var col = 0.99f * _gridSize * posGridSpace.x / GridWidth;

        Debug.Assert(row < _gridSize && col < _gridSize, $"Point {position} is outside of grid extents");

        return ((uint) Mathf.Floor(row), (uint) Mathf.Floor(col));
    }

    /// <summary>
    /// Return index of bin this point belongs to
    /// </summary>
    /// <param name="position"> Position in 2D space </param>
    /// <returns></returns>
    public uint BinIndex(Vector2 position)
    {
        var (row, col) = ComputeGridPosition(position);
        var result = row % 2 == 0 ? row * _gridSize + col + 1 : (row + 1) * _gridSize - col;

        return result -1;
    }

    public IEnumerator<Vector3> GetEnumerator()
    {
        foreach (var bin in _bins)
        {
            foreach (var vec in bin)
                yield return vec;
        }
    }

    IEnumerator IEnumerable.GetEnumerator()
    {
        throw new NotImplementedException();
    }
}
