using System.Collections.Generic;
using System.Collections;
using UnityEngine;

[RequireComponent(typeof(MeshFilter))]
public class procTreeScript : MonoBehaviour
{
    [Range(10, 5000)]
    public int numPoints;
    public int iteration;
    public float length;
    public bool orbs;


    Mesh mesh;
    public GameObject point;
    public GameObject volume;
    public LineRenderer line;
    public int number_of_points;

   
void Start()
    {
        //Physics.SyncTransforms();                       
        //GameObject volumeReal = (GameObject)Instantiate(volume, new Vector3(0.0f, dist, 0.0f), Quaternion.identity);
        //MeshCollider volumeCollider = volumeReal.GetComponent<MeshCollider>();
        //if (volumeCollider.bounds.Contains(position))
        
        TreeStruct test = new TreeStruct(0.5f, 4.0f, numPoints, length, number_of_points);

        test.growTree(iteration);

        //Debug.Log("number of leaves left: "+test.leafList.Count);
        //Debug.Log("total number of segments: "+test.segmentsList.Count);

        //test.draw(line,point);
        //test.drawLeaf(point);
        //test.calcVerts();
        //test.drawNewVerts(point);
        //test.drawVerts(point);
        //test.calcTris();

        test.createVerts();
        test.createMesh();
        test.CalculateNormals();
     
       

        mesh = new Mesh();
        GetComponent<MeshFilter>().mesh = mesh;
        mesh.vertices = test.verts;
        mesh.triangles = test.tris;
        mesh.normals = test.normals;

        mesh.uv = test.uvs;
        
        //Debug.Log(test.verts.Length);
        Debug.Log(test.tris.Length/3);

    }

    //Simple node like class, stores a position and a bool
    class Leaf
    {
        public Vector3 position;       
        public bool visted;
        public Leaf(Vector3 pos)
        {
            position = pos;
            visted = false;

        }       
    }

    //Class to store the information about each segment
    class Segment
    {
        public Vector3 start;
        public Vector3 end;
        public Vector3 direction;
        private Vector3 saveDir;
        public List<Vector3> startVerts;
        public List<Vector3> endVerts;
        public int count;
        public int age;
        public int index;

        public float segmentLength;

        public Segment parent;
        public List<Segment> children;


        public Segment(Vector3 srt, Vector3 dir,float len)
        {
            segmentLength = len;
            direction = dir;
            saveDir = dir;

            start = srt;
            end = start + (dir * segmentLength);
            
            count = 0;

            parent = null;
            children = new List<Segment>();
            age = 0;
        }

        public Segment(Segment p)
        {
            parent = p;
            segmentLength = p.segmentLength;
            direction = p.direction;
            saveDir = direction;

            start = p.end;
            end = start + (direction * segmentLength);

            count = 0;
            children = new List<Segment>();
            p.addChild(this);
            age = p.age + 1;
        }

        public void addChild(Segment c)
        {
            children.Add(c);
        }

        public Vector3 avgDir()
        {
            Vector3 output = direction;
            foreach (Segment seg in children)
            {
                output += seg.direction;
            }

            //output /= children.Count + 1;

            return output.normalized;
        }

        public void reset()
        {
            count = 0;
            direction = saveDir;

        }
    }

    //Class to store the leaves and segments, and to build the tree
    class TreeStruct
    {
        public List<Leaf> leafList;
        public List<Segment> segmentsList;
        public Vector3[] verts;
        public int[] tris;
        public Vector2[] uvs;
        public Vector3[] normals;
        private Dictionary<Vector3,List<int>> normalIndex;

        float maxDist;
        float minDist;
        int numLeaf;
        Segment root;
        int points;

        public TreeStruct(float min,float max,int nl, float length,int numPoints)
        {
            //Initalize variables
            maxDist = max;
            minDist = min;
            numLeaf = nl;
            points = numPoints;

            leafList = new List<Leaf>();
            segmentsList = new List<Segment>();
            normalIndex = new Dictionary<Vector3, List<int>>();

            float dist = 5;
            Vector3 center = new Vector3(0, 6, 0);


            //Loop until there are the correct number of leaves placed inside the volume
            while (leafList.Count < numLeaf)
            {
                Vector3 position = new Vector3(UnityEngine.Random.Range(-dist, dist), UnityEngine.Random.Range(-dist + 6, dist + 6), UnityEngine.Random.Range(-dist, dist));

                //Check if the new point is in the volume
                if (Vector3.Distance(position, center) < 3)
                {
                    leafList.Add(new Leaf(position));
                }
            }

            //Create the root node and its first child segment
            root = new Segment(new Vector3(0, 0, 0), new Vector3(0, 1, 0), length);
            segmentsList.Add(root);
            Segment current = new Segment(root);
            segmentsList.Add(current);

            //Grow the root until it reaches the clound of leaves
            while (!reached(current))
            {
                Segment seg = new Segment(current);
                segmentsList.Add(seg);
                current = seg;
            }
        }
        
        //Check if the root segments have reached the leaf cloud
        public bool reached(Segment seg)
        {
            foreach (Leaf l in leafList)
            {
                float dist = Vector3.Distance(seg.end, l.position);
                if (dist < maxDist)
                {
                    return true;
                }
            }
            return false;
        }

        public int updateTree()
        {
            //Loop through all leaves left in the tree
            foreach (Leaf l in leafList)
            {
                //Initalize starting values
                Segment closest = null;
                Vector3 closestDir = Vector3.zero;
                float currentDist = -1;

                //Loop through all segments in the tree
                foreach (Segment s in segmentsList)
                {
                    //calculate the length of the current leaf and segment
                    Vector3 dir = l.position - s.end;
                    float len = dir.magnitude;
                    
                    //Check if the length is less than the min distance
                    if (len < minDist)
                    {
                        //If less than min distance set the leaf to be removed and break from the loop
                        l.visted = true;
                        closest = null;
                        break;
                    
                    } else if(closest==null || len < currentDist)
                    {                       
                        closest = s;
                        closestDir = dir;
                        currentDist = len;
                    }
                }

                if (closest != null)
                {
                    closestDir.Normalize();
                    closest.direction += closestDir;
                    closest.count++;
                }
            }

            for (int i = leafList.Count-1; i > 0; i--)
            {                
                if (leafList[i].visted)
                {                    
                    leafList.RemoveAt(i);
                }
            }


            for (int i = segmentsList.Count-1; i > 0; i--)
            {
                Segment s = segmentsList[i];
                if (s.count>0)
                {
                    s.direction /= s.count;
                    s.direction.Normalize();

                    // Code top check if the next branch to be added is in the same place as an existing branch
                    bool childTest = true;
                    foreach (Segment child in s.children) { 
                        Vector3 tempVec = s.end + (s.direction * s.segmentLength);
                        if (child.end == tempVec)
                        {
                            childTest = false;
                            break;
                        }    
                    }

                    if (childTest)// && s.children.Count<2)
                    {
                        Segment seg = new Segment(s);
                        segmentsList.Add(seg);
                    }
                    s.reset();
                }
            }           
        
            return leafList.Count;         
        }

        public void growTree(int iter) {
            for (int i = 0; i < iter; i++) {
               updateTree();
            }
            cleanTree();
        }

        private void cleanTree() {
            for (int i = segmentsList.Count-1; i > 2; i--){
                Segment current = segmentsList[i];
                if (current.children.Count == 0) {
                    if (current.parent.parent.children.Count>1) {
                        current.parent.children.Remove(current);
                        segmentsList.RemoveAt(i);
                    }else if (current.parent.children.Count > 1) {
                        current.parent.children.Remove(current);
                        segmentsList.RemoveAt(i);
                    }
                }
            }
        }      

        public void createVerts() {

            List<Vector3> tempVerts = new List<Vector3>();
            
            float l = (4 / ((float)segmentsList[0].age + 4)) / 1.5f;

            tempVerts.Add(new Vector3(l, 0, 0));

            float numPoints = points;

            Vector3 axis = new Vector3(0, 360 / numPoints, 0);

            for (int i = 0; i < points - 1; i++) {
                tempVerts.Add(RotatePointAroundPivot(tempVerts[tempVerts.Count - 1], Vector3.zero, axis));
            }
            
            root.startVerts = tempVerts;           
           
            addToMesh(root);
        }

        private void addToMesh(Segment currentSeg) {
            //Code to add the verts and tris for current segment if cylinder
            if (currentSeg.children.Count >0) {
                getRingVerts(currentSeg);
                foreach (Segment child in currentSeg.children) {
                    
                    addToMesh(child);
                }

            //Code to add the verts and tris for current segment if end of a branch
            } else {
                currentSeg.startVerts = currentSeg.parent.endVerts;
                currentSeg.endVerts = new List<Vector3>();
                currentSeg.endVerts.Add(currentSeg.start+currentSeg.direction*currentSeg.segmentLength/4);
            }     
        }

        private void getRingVerts(Segment currentSeg) {

            if (currentSeg.parent != null) {
                currentSeg.startVerts = currentSeg.parent.endVerts;
            }

            List<Vector3> tempVerts = new List<Vector3>();

            float length = (4 / ((float)currentSeg.age + 4)) / 1.5f;
            float oldLength = (4 / ((float)currentSeg.age - 1 + 4)) / 1.5f;

            for (int i = 0; i < currentSeg.startVerts.Count; i++) {
                //create second ring at end of segment
                tempVerts.Add(currentSeg.startVerts[i] + (currentSeg.direction * currentSeg.segmentLength));
                //create the rotation matrix
                if (currentSeg.parent != null) {
                    Quaternion rotation = Quaternion.FromToRotation(currentSeg.parent.avgDir(), currentSeg.avgDir());

                    //rotate ring
                    tempVerts[i] = RotatePointAroundPivot(tempVerts[i], currentSeg.end, rotation);
                    //tempVerts[i] = (rotation * tempVerts[i]);
                }


                //create the temp vec for scaling
                Vector3 tempScale = tempVerts[i] - currentSeg.end;
                //scale vector
                float scale = length / oldLength;
                tempScale.Scale(new Vector3(scale, scale, scale));
                tempVerts[i] = tempScale + currentSeg.end;
            }

            currentSeg.endVerts = tempVerts;
        }

        public void createMesh() {
            List<Vector3>   tempVerts   = new List<Vector3>();
            List<int>       tempTris    = new List<int>();
            List<Vector2>   tempUV      = new List<Vector2>();

            int currentIndex = 0;
            float length = (float)1 / (float)points;
            Debug.Log(length);
            foreach (Segment currentSeg in segmentsList) {
                
                if (currentSeg.children.Count > 0) {

                    tempVerts.Add(currentSeg.startVerts[0]);
                    addToHashmap(currentSeg.startVerts[0], currentIndex);
                    tempTris.Add(currentIndex++);

                    tempVerts.Add(currentSeg.endVerts[0]);
                    addToHashmap(currentSeg.endVerts[0], currentIndex);
                    tempTris.Add(currentIndex++);

                    tempVerts.Add(currentSeg.endVerts[0+points-1]);
                    addToHashmap(currentSeg.endVerts[0 + points - 1], currentIndex);
                    tempTris.Add(currentIndex++);


                    tempUV.Add(new Vector2(1, 0));
                    tempUV.Add(new Vector2(1, 1));
                    tempUV.Add(new Vector2(1-length, 1));


                    tempVerts.Add(currentSeg.startVerts[0]);
                    addToHashmap(currentSeg.startVerts[0], currentIndex);
                    tempTris.Add(currentIndex++);

                    tempVerts.Add(currentSeg.endVerts[0 + points - 1]);
                    addToHashmap(currentSeg.endVerts[0 + points - 1], currentIndex);
                    tempTris.Add(currentIndex++);

                    tempVerts.Add(currentSeg.startVerts[0 + points - 1]);
                    addToHashmap(currentSeg.startVerts[0 + points - 1], currentIndex);
                    tempTris.Add(currentIndex++);

                    tempUV.Add(new Vector2(1, 0));
                    tempUV.Add(new Vector2(1-length, 1));
                    tempUV.Add(new Vector2(1-length, 0));

                    for (int j = 1; j < points; j++) {
                        tempVerts.Add(currentSeg.startVerts[0+j]);
                        addToHashmap(currentSeg.startVerts[0+j], currentIndex);
                        tempTris.Add(currentIndex++);

                        tempVerts.Add(currentSeg.endVerts[0+j]);
                        addToHashmap(currentSeg.endVerts[0+j], currentIndex);
                        tempTris.Add(currentIndex++);

                        tempVerts.Add(currentSeg.endVerts[0 + j - 1]);
                        addToHashmap(currentSeg.endVerts[0 + j - 1], currentIndex);
                        tempTris.Add(currentIndex++);

                        tempUV.Add(new Vector2(length * j, 0));
                        tempUV.Add(new Vector2(length * j, 1));
                        tempUV.Add(new Vector2(length * (j - 1), 1));


                        tempVerts.Add(currentSeg.startVerts[0 + j]);
                        addToHashmap(currentSeg.startVerts[0 + j], currentIndex);
                        tempTris.Add(currentIndex++);

                        tempVerts.Add(currentSeg.endVerts[0 + j - 1]);
                        addToHashmap(currentSeg.endVerts[0 + j - 1], currentIndex);
                        tempTris.Add(currentIndex++);

                        tempVerts.Add(currentSeg.startVerts[0 + j - 1]);
                        addToHashmap(currentSeg.startVerts[0 + j - 1], currentIndex);
                        tempTris.Add(currentIndex++);

                        tempUV.Add(new Vector2(length * j, 0));
                        tempUV.Add(new Vector2(length * (j-1), 1));
                        tempUV.Add(new Vector2(length * (j-1), 0));
                    }

                } else {
                    tempVerts.Add(currentSeg.startVerts[0]);
                    addToHashmap(currentSeg.startVerts[0], currentIndex);
                    tempTris.Add(currentIndex++);

                    tempVerts.Add(currentSeg.endVerts[0]);
                    addToHashmap(currentSeg.endVerts[0], currentIndex);
                    tempTris.Add(currentIndex++);

                    tempVerts.Add(currentSeg.startVerts[0 + points - 1]);
                    addToHashmap(currentSeg.startVerts[0 + points - 1], currentIndex);
                    tempTris.Add(currentIndex++);

                    tempUV.Add(new Vector2(1, 0));
                    tempUV.Add(new Vector2(1 - ((length / 2)), 1));
                    tempUV.Add(new Vector2(1 - length, 0));

                    //tempUV.Add(new Vector2(0, 0));
                    //tempUV.Add(new Vector2(0.5f, 1));
                    //tempUV.Add(new Vector2(1, 0));

                    for (int j = 1; j < points; j++) {
                        tempVerts.Add(currentSeg.startVerts[0+ j]);
                        addToHashmap(currentSeg.startVerts[0 + j], currentIndex);
                        tempTris.Add(currentIndex++);

                        tempVerts.Add(currentSeg.endVerts[0]);
                        addToHashmap(currentSeg.endVerts[0], currentIndex);
                        tempTris.Add(currentIndex++);

                        tempVerts.Add(currentSeg.startVerts[0 + j - 1]);
                        addToHashmap(currentSeg.startVerts[0 + j - 1], currentIndex);
                        tempTris.Add(currentIndex++);

                        tempUV.Add(new Vector2(length * j, 0));
                        tempUV.Add(new Vector2((length * (j - 1)) + length / 2, 1));
                        tempUV.Add(new Vector2(length * (j - 1), 0));

                        //tempUV.Add(new Vector2(0, 0));
                        //tempUV.Add(new Vector2(0.5f, 1));
                        //tempUV.Add(new Vector2(1, 0));
                    }
                }
            }
            
            verts = tempVerts.ToArray();
            tris = tempTris.ToArray();
            uvs = tempUV.ToArray();
        }

        private void addToHashmap(Vector3 vec, int index) {
            if (!normalIndex.ContainsKey(vec)) {
                normalIndex.Add(vec, new List<int>());
                normalIndex[vec].Add(index);
            } else {
                normalIndex[vec].Add(index);
            }
        }

        public void CalculateNormals() {
            normals = new Vector3[tris.Length];
            for (int i = 0; i < tris.Length; i += 3) {
                int tri0 = tris[i];
                int tri1 = tris[i + 1];
                int tri2 = tris[i + 2];
                Vector3 vert0 = verts[tri0];
                Vector3 vert1 = verts[tri1];
                Vector3 vert2 = verts[tri2];
                // Vector3 normal = Vector3.Cross(vert1 - vert0, vert2 - vert0);
                Vector3 normal = new Vector3() {
                    x = vert0.y * vert1.z - vert0.y * vert2.z - vert1.y * vert0.z + vert1.y * vert2.z + vert2.y * vert0.z - vert2.y * vert1.z,
                    y = -vert0.x * vert1.z + vert0.x * vert2.z + vert1.x * vert0.z - vert1.x * vert2.z - vert2.x * vert0.z + vert2.x * vert1.z,
                    z = vert0.x * vert1.y - vert0.x * vert2.y - vert1.x * vert0.y + vert1.x * vert2.y + vert2.x * vert0.y - vert2.x * vert1.y
                };
                normals[tri0] += normal;
                normals[tri1] += normal;
                normals[tri2] += normal;
            }

            for (int i = 0; i < normals.Length; i++) {
                // normals [i] = Vector3.Normalize (normals [i]);
                Vector3 norm = normals[i];
                float invlength = 1.0f / (float)System.Math.Sqrt(norm.x * norm.x + norm.y * norm.y + norm.z * norm.z);
                normals[i].x = norm.x * invlength;
                normals[i].y = norm.y * invlength;
                normals[i].z = norm.z * invlength;
            }

            foreach (KeyValuePair<Vector3, List<int>> list in normalIndex) {
                Vector3 newNormal = Vector3.zero;
                foreach (int index in list.Value) {
                    newNormal += normals[index];
                }
                newNormal.Normalize();
                foreach (int index in list.Value) {
                    normals[index] = newNormal;
                }
            }
        }

        public void calcVerts()
        {
            int currentindex = 5;
            List<Vector3> tempVerts = new List<Vector3>();

            float l = (4 / ((float)segmentsList[0].age + 4)) / 1.5f;

            tempVerts.Add(new Vector3(l, 0, 0));            

            float points = 5;

            Vector3 axis = new Vector3(0, 360 / points, 0);

            for (int i = 0; i < points - 1; i++)
            {
                tempVerts.Add(RotatePointAroundPivot(tempVerts[tempVerts.Count - 1], Vector3.zero, axis));
            }

            for (int i = tempVerts.Count - ((int)points); i < tempVerts.Count; i++)
            {
                Quaternion rotation = Quaternion.FromToRotation(Vector3.up, segmentsList[0].avgDir());
                tempVerts[i] = (rotation * tempVerts[i]) + segmentsList[0].start;
            }

            segmentsList[0].index = currentindex;


            for (int j = 0; j < segmentsList.Count; j++)
            {
                Segment s = segmentsList[j];
                if (s.children.Count == 0)
                {
                    tempVerts.Add(s.end);
                    
                    //might need to check for if j is 0
                    if (segmentsList[j-1].children.Count==0) {
                        currentindex += 1;
                        s.index = currentindex;
                    } else {
                        currentindex += 5;
                        s.index = currentindex;
                    }
                    
                }
                else 
                {

                    l = (4 / ((float)s.age + 4)) / 1.5f;

                    tempVerts.Add(new Vector3(l, 0, 0));                                      

                    for (int i = 0; i < points - 1; i++)
                    {
                        tempVerts.Add(RotatePointAroundPivot(tempVerts[tempVerts.Count - 1], Vector3.zero, axis));
                    }

                    for (int i = tempVerts.Count - ((int)points); i < tempVerts.Count; i++)
                    {
                        Quaternion rotation = Quaternion.FromToRotation(Vector3.up, s.avgDir());
                        tempVerts[i] = (rotation * tempVerts[i]) + s.end;
                    }

                    if (j != 0) {
                        if (segmentsList[j - 1].children.Count == 0) {
                            currentindex += 1;
                            s.index = currentindex;
                        } else {
                            currentindex += 5;
                            s.index = currentindex;
                        }
                    }
                }
                
            }
            verts = tempVerts.ToArray();

            Debug.Log("total verts: " + verts.Length);
            Debug.Log("last index: " + currentindex);
        }

        public void calcTris() {
            List<int> tempTris = new List<int>();
            List<Vector2> tempUV = new List<Vector2>();

            for (int i = 1; i<segmentsList.Count-1; i++) {
                for (int child = 0; child < segmentsList[i].children.Count; child++) {
                    if (segmentsList[i].children[child].children.Count != 0) {
                        tempTris.Add(segmentsList[i].index);
                        tempTris.Add(segmentsList[i].children[child].index);
                        tempTris.Add(segmentsList[i].children[child].index + 4);                            

                        tempTris.Add(segmentsList[i].index);                        
                        tempTris.Add(segmentsList[i].children[child].index + 4);
                        tempTris.Add(segmentsList[i].index + 4);

                        for (int j = 1; j < 5; j++) {
                            tempTris.Add(segmentsList[i].index+j);
                            tempTris.Add(segmentsList[i].children[child].index + j);
                            tempTris.Add(segmentsList[i].children[child].index + j - 1);

                            tempTris.Add(segmentsList[i].index + j);
                            tempTris.Add(segmentsList[i].children[child].index + j - 1);
                            tempTris.Add(segmentsList[i].index + j - 1);
                        }
                    } else {

                        tempTris.Add(segmentsList[i].index);
                        tempTris.Add(segmentsList[i].children[child].index);
                        tempTris.Add(segmentsList[i].index + 4);

                        tempUV.Add(new Vector2(0, 0));
                        tempUV.Add(new Vector2(0.5f, 1));
                        tempUV.Add(new Vector2(1, 0));

                        for (int j = 1; j < 5; j++) {
                            tempTris.Add(segmentsList[i].index + j);
                            tempTris.Add(segmentsList[i].children[child].index);
                            tempTris.Add(segmentsList[i].index + j - 1);

                            tempUV.Add(new Vector2(0, 0));
                            tempUV.Add(new Vector2(0.5f, 1));
                            tempUV.Add(new Vector2(1, 0));
                        }
                    }
                }
            }

            tris = tempTris.ToArray();
            uvs = tempUV.ToArray();
        }

        public Vector3 RotatePointAroundPivot(Vector3 point, Vector3 pivot, Vector3 angles)
        {
            Vector3 dir = point - pivot;            // get point direction relative to pivot
            dir = Quaternion.Euler(angles) * dir;   // rotate it
            point = dir + pivot;                    // calculate rotated point
            return point;                           
        }

        public Vector3 RotatePointAroundPivot(Vector3 point, Vector3 pivot, Quaternion matrix) {
            Vector3 dir = point - pivot;            // get point direction relative to pivot
            dir = matrix * dir;   // rotate it
            point = dir + pivot;                    // calculate rotated point
            return point;
        }

        public void draw(LineRenderer line, GameObject point)
        {
            foreach (Segment s in segmentsList)
            {
                LineRenderer seg = Instantiate(line, new Vector3(0, 0, 0), Quaternion.identity);
                seg.SetPosition(0, s.start);
                seg.SetPosition(1, s.end);              
            }
        }

        public void drawLeaf(GameObject point)
        {
            foreach (Leaf s in leafList)
            {
                GameObject temp = Instantiate(point, s.position, Quaternion.identity);              
            }
        }

        public void drawVerts(GameObject point)
        {
            foreach (Vector3 vert in verts)
            {
                GameObject temp = Instantiate(point, vert, Quaternion.identity);
            }
        }

        public void drawNewVerts(GameObject point) { 
            foreach(Segment current in segmentsList) {
                foreach(Vector3 svert in current.startVerts) {
                    GameObject temp = Instantiate(point, svert, Quaternion.identity);
                }

                foreach (Vector3 evert in current.endVerts) {
                    GameObject temp = Instantiate(point, evert, Quaternion.identity);
                }
            }
        }
    }
}