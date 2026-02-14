function tri=makeTriangle(w,h,yLoc)
    x=[0 w w];%x coordinates of vertices
    y=[yLoc -h+yLoc h+yLoc];%y coordinates of vertices
    tri = patch(x,y,'black');
end