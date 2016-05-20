r.proj location=colimafine mapset=PERMANENT input=colima
r.fillnulls --overwrite input=colima@PERMANENT output=colima method=rst
