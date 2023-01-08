using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[RequireComponent(typeof(MeshFilter))]
public class DeluanayTriangulationTestComponent : MonoBehaviour
{
    // Mesh to start shape from
    private MeshFilter _meshFilter;

    // List of points added by the user
    private List<Vector2> _newPoints;

    // List of points corresponding to holes
    private List<Vector2> _holePoints;

    private PointBinGrid _pointBinGrid;

    private Dictionary<uint, Color> _binIndexToColor;

    private DeluanayTriangulation.TriangulationResult? _result;

    // Start is called before the first frame update
    void Start()
    {
        _meshFilter = GetComponent<MeshFilter>();
        Debug.Assert(_meshFilter != null, "Mesh component missing");
        Debug.Assert(_meshFilter.mesh.vertices.Length == 4, "Mesh should be a quad");
        _newPoints = new List<Vector2>();
        _holePoints = new List<Vector2>();
        _binIndexToColor = new Dictionary<uint, Color>();
    }

    // Update is called once per frame
    void Update()
    {
        if (Input.GetMouseButtonDown(0))
            TryAddNewPoint();
        else if (Input.GetMouseButtonDown(1))
            TryAddNewPoint(true);
        else if (Input.GetKeyDown(KeyCode.Space))
        {
            _result = DeluanayTriangulation.Triangulate(_newPoints.ToArray(), new List<Vector2[]>{_holePoints.ToArray()});
            Debug.Log("Issuing triangulation request...");
        }
        else if (Input.GetKeyDown(KeyCode.C))
        {
            _newPoints.Clear();
            _holePoints.Clear();
        }    
    }


    void OnDrawGizmos()
    {
        if (_newPoints == null || _newPoints.Count == 0)
            return;

        Gizmos.color = Color.red;
        foreach (var point in _newPoints)
        {
            Gizmos.DrawSphere(point, 0.1f);
        }

        Gizmos.color = Color.blue;
        foreach(var point in _holePoints)
        {
            Gizmos.DrawSphere(point, 0.1f);
        }

        if (_result == null || _holePoints.Count > 0)
            return;

        Gizmos.color = Color.red;
        var actualResult = (DeluanayTriangulation.TriangulationResult)_result;
        foreach(var point in actualResult.points)
        {
            Gizmos.DrawSphere(point, 0.1f);
        }

        for(int i = 0; i < actualResult.triangles.Length; i += 3)
        {
            var p1 = actualResult.points[actualResult.triangles[i]];
            var p2 = actualResult.points[actualResult.triangles[i+1]];
            var p3 = actualResult.points[actualResult.triangles[i+2]];

            Gizmos.DrawLine(p1, p2);
            Gizmos.DrawLine(p2, p3);
            Gizmos.DrawLine(p3, p1);
        }
    }

    /// <summary>
    /// Try to add a new point based on mouse input if inside the right region
    /// </summary>
    private void TryAddNewPoint(bool isHole = false)
    {
        var maybeNewPoint = GetNewPointFromMouse();

        // If no hit, just exit
        if (maybeNewPoint == null)
            return;

        Vector3 newPoint = maybeNewPoint ?? Vector3.zero;

        // Check if it's  hole before adding it to the array of points
        if (isHole)
            _holePoints.Add(newPoint);
        else
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
