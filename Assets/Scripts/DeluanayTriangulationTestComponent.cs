using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[RequireComponent(typeof(MeshFilter))]
public class DeluanayTriangulationTestComponent : MonoBehaviour
{
    // Mesh to start shape from
    private MeshFilter _meshFilter;

    // List of points added by the user
    private List<Vector3> _newPoints;

    private PointBinGrid _pointBinGrid;

    private Dictionary<uint, Color> _binIndexToColor;

    // Start is called before the first frame update
    void Start()
    {
        _meshFilter = GetComponent<MeshFilter>();
        Debug.Assert(_meshFilter != null, "Mesh component missing");
        Debug.Assert(_meshFilter.mesh.vertices.Length == 4, "Mesh should be a quad");
        _newPoints = new List<Vector3>();
        _binIndexToColor = new Dictionary<uint, Color>();
    }

    // Update is called once per frame
    void Update()
    {
        if (Input.GetMouseButtonDown(0))
            TryAddNewPoint();
        else if (Input.GetKeyDown(KeyCode.Space))
            CreatePointBin();
    }

    void CreatePointBin()
    {
        if (_newPoints.Count <= 2)
        {
            Debug.LogWarning($"Can't create point bin, too few elements: {_newPoints.Count}");
            return;
        }

        // Assume that all points lie in the same plane
        var world2Local = transform.worldToLocalMatrix;
        Vector2[] points2d = new Vector2[_newPoints.Count];

        uint i = 0;
        foreach(var point in _newPoints)
        {
            var pointHomoCoords = new Vector4(point.x, point.y, point.z, 1.0f);
            var pointLocalHomoCoords = world2Local * pointHomoCoords;

            points2d[i++] = new Vector2(pointLocalHomoCoords.x, pointLocalHomoCoords.y);
        }

        // Create point bin from current points
        _pointBinGrid = new PointBinGrid(points2d);
    }

    void OnDrawGizmos()
    {
        if (_newPoints == null)
            return;

        Gizmos.color = Color.red;
        foreach (var point in _newPoints)
        {
            Gizmos.DrawSphere(point, 0.1f);
        }

        if (_pointBinGrid == null)
            return;

        uint i = 0;
        foreach (var bin in _pointBinGrid.Bins)
        {
            Color color;
            if (!_binIndexToColor.TryGetValue(i, out color))
            {
                color = new Color(
                    Random.Range(0.0f,1.0f),
                    Random.Range(0.0f,1.0f),
                    Random.Range(0.0f,1.0f)
                );

                _binIndexToColor[i] = color;
            }
            i++;  
            Gizmos.color = color;

            foreach(var point in bin)
            {
                // We have to convert from 2d point to 3d and then to world coordinates
                var pointHomoLocalCoords = new Vector4(point.x, point.y, 0.0f, 1.0f);
                var pointHomoWorldCoords = transform.localToWorldMatrix * pointHomoLocalCoords;
                var pointWorld3D = new Vector3(pointHomoWorldCoords.x, pointHomoWorldCoords.y, pointHomoWorldCoords.z);

                Gizmos.DrawSphere(pointWorld3D, 0.1f);
            }

        }
    }

    /// <summary>
    /// Try to add a new point based on mouse input if inside the right region
    /// </summary>
    private void TryAddNewPoint()
    {
        var maybeNewPoint = GetNewPointFromMouse();

        // If no hit, just exit
        if (maybeNewPoint == null)
            return;

        Vector3 newPoint = maybeNewPoint ?? Vector3.zero;

        _newPoints.Add(newPoint);

        Debug.Log($"Addning new point. Current point count: {_newPoints.Count}");
    }

    /// <summary>
    /// Try to get a new point to add from the mouse
    /// </summary>
    /// <returns> Null if no point was added to the current extents </returns>
    private Vector3? GetNewPointFromMouse()
    {
        var ray = Camera.main.ScreenPointToRay(Input.mousePosition);
        RaycastHit hit;
        if (Physics.Raycast(ray, out hit) && hit.collider.gameObject == gameObject)        
            return hit.point;

        return null;
    }
}
