package ndm.NetworkAnalyse;

import java.sql.CallableStatement;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Struct;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.Hashtable;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import org.apache.log4j.Logger;



import ndm.domain.ClosestPoint;
import ndm.domain.ClosestPointAndDist;
import ndm.domain.Tra_Node;
import oracle.jdbc.OracleTypes;
import oracle.spatial.geometry.JGeometry;
import oracle.spatial.network.Network;
import oracle.spatial.network.NetworkManager;
import oracle.spatial.network.Path;
import oracle.spatial.network.SubPath;


public class Match {
	
	private static Logger logger = Logger.getLogger(Match.class);
	Tra_Node[] tra; //修改后的轨迹数据
	double bufferdist = 0; //缓冲距离
	double adjacentdist = 0; //临近距离
	double limitspeed = 0;  //限制速度
	double limittime = 0;  //两点之间的最大时间差
	Connection conn; //数据库连接类
	int constvalue = 100; //防止概率太小被忽略 
	int tra_fileid;
	public Match(Connection conn,int tra_fileid,int feacount,double bufferdist,double adjacentdist,double limitspeed,double limittime)
	{
		tra = new Tra_Node[feacount];
		this.bufferdist = bufferdist;
		this.adjacentdist = adjacentdist;
		this.limitspeed = limitspeed;
		this.limittime = limittime;
		this.conn = conn;	
		this.tra_fileid = tra_fileid;
	}
	
	/**
	 * @param tra_fileid 文件ID
	 */
	public boolean MapMatch(Network wh_roads)
	{
		try 
		{	
			//Network wh_roads = NetworkManager.readNetwork(conn, "wh_roads"); //加载网络
			List<List<Tra_Node>> pTree = new ArrayList<List<Tra_Node>>();
			String sql = "select id,fileid,geo,time from wh_tras where fileid = ? order by id";
			PreparedStatement ps = conn.prepareStatement(sql,ResultSet.TYPE_SCROLL_INSENSITIVE, ResultSet.CONCUR_UPDATABLE);
			ps.setInt(1, tra_fileid);
			ResultSet rs_tra = ps.executeQuery(); //获取当前轨迹要素集
			List<Tra_Node> invalidPoint = new ArrayList<Tra_Node>(); //用于存储轨迹点缓冲区范围内没有道路的点
			List<Tra_Node> timeExceptionPoint = new ArrayList<Tra_Node>(); //用于存储时间异常的轨迹点
			List<AdjacnetPoint> adjacentPoint = new ArrayList<AdjacnetPoint>();//用于存储轨迹点与前一点相差太近的点
			Hashtable<Integer, List<Tra_Node>> allNodes = new Hashtable<Integer, List<Tra_Node>>();
			while (rs_tra.next())
			{	
				List<Tra_Node> gpsNodes = new ArrayList<Tra_Node>(); //用于存储当前点可能的道路
				//获取当前轨迹点对象
	            Struct struct_point = (Struct)rs_tra.getObject ("geo");
	            JGeometry geom_point = JGeometry.loadJS ( struct_point );
	            //System.out.println(geom_point.getPoint()[0]);
	            int tra_id = rs_tra.getInt("id");
	            //System.out.println(tra_id);
	            int fileid = rs_tra.getInt("fileid");
	            String time = rs_tra.getString("time");
	            //获取缓冲区内道路
	            String sql2 = "select l1.link_id,l1.street_geom,l1.start_node_id,l1.end_node_id,l1.street_length,l1.BIDIRECTED from links l1,wh_tras tra where "
	            		+ "tra.id= ? and tra.fileid= ? and sdo_within_distance(l1.street_geom,"
	            		+ "tra.geo,'DISTANCE=50 UNIT=METER')='TRUE'";
	            PreparedStatement ps2 = conn.prepareStatement(sql2,ResultSet.TYPE_SCROLL_INSENSITIVE, ResultSet.CONCUR_UPDATABLE);
	            ps2.setInt(1, tra_id);
			    ps2.setInt(2, tra_fileid);
			    //ps.setDouble(3, bufferdist);
			    ResultSet rs_links = ps2.executeQuery();	
			    int flag = 0; //用于判断rs是否为空
			    List<Integer> tra_links = new ArrayList<Integer>(); //用于保存缓冲区内的道路
			    while(rs_links.next())
			    {
			    	flag++;
			    	int link_id = rs_links.getInt("link_id");
			    	tra_links.add(link_id);
			    	Struct struct_link = (Struct)rs_links.getObject ("street_geom");
			    	ClosestPoint cpt = GetCPTAndDist(struct_point, struct_link); //获取最近点和距离
			    	int start_node_id = rs_links.getInt("start_node_id");
			    	int end_node_id = rs_links.getInt("end_node_id");
			    	String oneway = rs_links.getString("BIDIRECTED");
			    	JGeometry geo_link = JGeometry.loadJS(struct_link);
			    	double link_length = rs_links.getDouble("street_length");
			    	double startper = GetPercentage("wh_roads2", link_id, cpt.struct);
			    	//观测概率
			    	double mp_cost = 1.0 / (Math.sqrt(2 * Math.PI) * 9) * Math.exp(-0.5 * Math.pow((cpt.dist / 9), 2)) * constvalue;
			    	Tra_Node pNode = new Tra_Node(-1, tra_id, link_id, fileid, time, 
			    			mp_cost,struct_point,geom_point,cpt.struct,cpt.geo,struct_link,geo_link,
			    			start_node_id,end_node_id,oneway,link_length,startper,null,null,-1);
			    	gpsNodes.add(pNode); //*****注意root_index、cost、path在下面需要修改
			    	//allNodes.put(tra_id, gpsNodes);
			    }
			    allNodes.put(tra_id, gpsNodes);
			    rs_links.close();
			    int _index = pTree.size();
//			    //*********************处理轨迹点突然偏移，正确道路不在范围内的情况***********************
//			    if (_index != 0)
//			    {
//			    	int beforelinkid = pTree.get(_index-1).get(0).link_id;
//			    	double cost_temp =  pTree.get(_index-1).get(0).cost;
//			    	int nodeid_temp = pTree.get(_index-1).get(0).link_endNodeID;
//			    	for (int i = 1; i < _index; i++)
//			    	{
//						if (pTree.get(_index-1).get(i).cost>cost_temp)
//						{
//							beforelinkid = pTree.get(_index-1).get(i).link_id;
//							cost_temp = pTree.get(_index-1).get(i).cost;
//							nodeid_temp = pTree.get(_index-1).get(i).link_endNodeID;
//						}
//					}
//					String sql_link = " select l1.link_id,l1.street_geom,l1.start_node_id,l1.end_node_id,"
//							+ "l1.street_length,l1.BIDIRECTED from links l1 where l1.link_id = ? or "
//							+ "l1.start_node_id = ? or l1.end_node_id = ?";
//					ps = conn.prepareStatement(sql_link,ResultSet.TYPE_SCROLL_INSENSITIVE, ResultSet.CONCUR_UPDATABLE);
//				    ps.setInt(1, beforelinkid);
//				    ps.setInt(2, nodeid_temp);
//				    ps.setInt(3, nodeid_temp);
//				    ResultSet rs_temp = ps.executeQuery();	
//				    while (rs_temp.next())
//				    {
//				    	int link_id = rs_temp.getInt("link_id");
//				    	if (tra_links.contains(link_id))
//				    	{
//							continue;
//						}
//				    	
//						
//					}
//			    }
//			    //*************************************************************
			    ps2.close();
	            if (flag == 0)
	            {
	            	invalidPoint.add(new Tra_Node(-1, tra_id, -1, fileid, time,
	            			-1, struct_point,geom_point, null, null, null, null, -1, -1, null, -1, -1, null,null,-1));
					continue; //如果为空，说明当前轨迹点缓冲区范围内没有道路，下一个  *****
				}	 
	            //int _index = pTree.size();
	            if (_index!=0) 
	            {
	            	Date date = GetDate(time);
	            	Date date_before = GetDate(pTree.get(_index-1).get(0).time);
	            	long timespan = Math.abs(date.getTime() - date_before.getTime())/1000;
	            	if (timespan > limittime*60)
	            	{
	            		//*****************************************
	            		WriteRoute(pTree, invalidPoint, adjacentPoint);
	            		pTree.clear();
	            		invalidPoint.clear();
	            		//timeExceptionPoint.clear();
	            		adjacentPoint.clear();
	            		rs_tra.previous();
	            		continue;//超过限定时间  此处应形成一个break					 	*****
					}	               
				    //double dist_Adjacent = geom_point.distance(pTree.get(_index-1).get(0).geo_old, 0.05, "true");
				    double dist_Adjacent = getDistance(struct_point, pTree.get(_index-1).get(0).st_old, conn);
				    //************************************
				    if (dist_Adjacent == 0)
				    {
						System.out.println("计算距离出现异常！135行");
						logger.debug("计算距离出现异常！135行");
					}
				    //************************************
				    
				    if (dist_Adjacent < adjacentdist)
					{
						adjacentPoint.add(new AdjacnetPoint(pTree.get(_index-1).get(0).id, struct_point,geom_point,tra_id, fileid, time) ); //此处只保存上一结点的ID，这样就可以通过modifiedGPS数组将结果移植过来
						continue; //小于设置的相邻距离      *****此处需要处理，将此点结果设置为上一点的
					}
					//将时间异常的点提取出来  **************** timespan < 5 是否正确？
					if ((dist_Adjacent/timespan) > (limitspeed/3.6) && timespan < 10) 
					{
						timeExceptionPoint.add(new Tra_Node(-1, tra_id, -1, fileid, time,
		            			-1, struct_point,geom_point, null, null, null, null, -1, -1, null, -1, -1, null,null,-1));
						continue;
					}
					List<Tra_Node> invalidgpsNode =new ArrayList<Tra_Node>();//将无效的node点的索引保存
					for(int j=0;j< gpsNodes.size(); j++ )
					{		
						double cost = 0;
						boolean invalidflag = false;
						for(int i = 0; i< pTree.get(_index-1).size();i++)
						{
							Tra_Node before = pTree.get(_index-1).get(i);
							SubPath subpath = null; //最短路径
							double dist_route = 0; //路径距离
							if (NetworkManager.isReachable(wh_roads, before.link_startNodeID, gpsNodes.get(j).link_endNodeID))
							{
								if (before.link_id == gpsNodes.get(j).link_id)
								{
									if (before.oneway.equals("N"))
									{
										if (gpsNodes.get(j).percentage-before.percentage<0)
										{
											dist_route = 0;
										}
										else
										{
											dist_route = (gpsNodes.get(j).percentage-before.percentage) * before.link_length;
										}
									}
									else
									{
										dist_route = Math.abs(gpsNodes.get(j).percentage-before.percentage) * before.link_length;
									}
								}
								else
								{
									//****************************************************
									if ((before.percentage>0 && before.percentage<1)||(gpsNodes.get(j).percentage>0 && gpsNodes.get(j).percentage<1))
									{
										subpath = GetSubPath(wh_roads, before, gpsNodes.get(j),timespan*(limitspeed/3.6));
									}
									else
									{
										int startNodeID = -1;
										if (before.percentage==0)
										{
											startNodeID = before.link_startNodeID;
										}
										else
										{
											startNodeID = before.link_endNodeID;
										}
										int endNodeID = -1;
										if (gpsNodes.get(j).percentage == 0)
										{
											endNodeID = gpsNodes.get(j).link_startNodeID;
										}
										else 
										{
											endNodeID = gpsNodes.get(j).link_endNodeID;
										}
										if (startNodeID != -1 && endNodeID != -1) 
										{
											if (startNodeID == endNodeID)
											{
												dist_route = 0;
											}
											else 
											{
												Path path = NetworkManager.shortestPath(wh_roads, startNodeID, endNodeID);
												if (path != null)
												{
													dist_route = path.getCost();
												}
											}										
										}
									}
									
								}
							}
							//求路径距离
							if (subpath != null)
							{
								dist_route = subpath.getCost();
							}
							//路径距离为0 或者超速 下一个
							if (dist_route == 0 || (dist_route/timespan)>(limitspeed/3.6))
							{
								continue;
							}
							invalidflag = true; //此点有路可达 有效
							double dt = Math.abs(dist_Adjacent-dist_route);
							double tp_cost = 0.039 * Math.exp(-dt * 0.039) * constvalue; //转移概率
							if (cost==0)
							{
								cost= gpsNodes.get(j).cost * tp_cost* before.cost;
								//gpsNodes.get(j).cost = cost;
								gpsNodes.get(j).root_index = i;
								if (subpath != null)
								{
									gpsNodes.get(j).path = subpath.getReferencePath();
								}
								else
								{
									gpsNodes.get(j).path = null;
								}							
								gpsNodes.get(j).subPath = subpath;
								gpsNodes.get(j).beforeid = before.id;
							}
							else
							{
								double temp = gpsNodes.get(j).cost * tp_cost* before.cost;
								if (temp > cost)
								{
									cost = temp;
									//gpsNodes.get(j).cost = cost;
									gpsNodes.get(j).root_index = i;
									if (subpath != null)
									{
										gpsNodes.get(j).path = subpath.getReferencePath();
									}
									else
									{
										gpsNodes.get(j).path = null;
									}
									gpsNodes.get(j).subPath = subpath;
									gpsNodes.get(j).beforeid = before.id;
								}
							}
						}
						gpsNodes.get(j).cost = cost;
						//添加无效点的索引
						if (!invalidflag)
						{
							invalidgpsNode.add(gpsNodes.get(j));
						}
						//System.out.println(gpsNodes.get(j).root_index);
					}
					//删除无效的GPS点
					if (invalidgpsNode.size()!=0)
					{
						gpsNodes.removeAll(invalidgpsNode);
					}
				}
	            if (gpsNodes.size()!=0)
	            {
					pTree.add(gpsNodes);
				}  
	            else //若个数为0，说明存在前后两点无法到达
	            {
	            	//*******************************************
	            	WriteRoute(pTree, invalidPoint, adjacentPoint);
            		pTree.clear();
            		invalidPoint.clear();
            		//timeExceptionPoint.clear();
            		adjacentPoint.clear();
	            	rs_tra.previous();
	            	continue;
	            }
			}
			rs_tra.close();
			ps.close();
			WriteRoute(pTree, invalidPoint, adjacentPoint);
			SolveTimeExceptionPoint(allNodes,timeExceptionPoint,wh_roads);
			return true;
		} 
		catch (Exception e)
		{
			e.printStackTrace(); 
			return false;
		}
	}
	
	/**
	 * @param wh_roads 网络数据集
	 * @param before   前一结点
	 * @param current  当前结点
	 * @param MaxCost  最大的距离
	 * @return 最短路径
	 */
	private SubPath GetSubPath(Network wh_roads,Tra_Node before,Tra_Node current,double MaxCost)
	{
		SubPath subpath = null;
		try 
		{
			double startPercentage = before.percentage;
			double endPercentage = current.percentage;
			subpath = NetworkManager.shortestPath(wh_roads, before.link_id,startPercentage,current.link_id, 
					endPercentage, null);		//new CostConstraint(MaxCost)	
			return subpath;			
		} 
		catch (Exception e) //NetworkDataException
		{
			//e.printStackTrace();
			System.out.println(e.getMessage());
			logger.error(e.getMessage());
			return null;
		}
	}
	
	/**
	 * @param pTree 所有结点
	 * @param invalidPoint 缓冲区内没有道路的点
	 * @param timeExceptionPoint 时间异常的点
	 * @param adjacentPoint 相邻两点距离太短的点
	 */
	private void WriteRoute(List<List<Tra_Node>> pTree,List<Tra_Node>invalidPoint,List<AdjacnetPoint> adjacentPoint) throws SQLException
	{		
		if (pTree.size()==0)
		{
			return;
		}
		//获取最后一点cost值最大的和索引
		 int pNode_Index = 0;
		 double cost = pTree.get(pTree.size()-1).get(0).cost;
		 for(int i=1;i<pTree.get(pTree.size()-1).size();i++)
		 {
			 if (pTree.get(pTree.size()-1).get(i).cost>cost)
			 {
				pNode_Index = i;
				cost = pTree.get(pTree.size()-1).get(i).cost;
			 }
		 }
		 
		 String sql_new = "insert into wh_modifiedTra(new,tra_id,file_id,time,roadid) values(?,?,?,?,?)"; 
		 String sql_old = "insert into wh_tras_old(origin,tra_id,file_id,time,road_id) values(?,?,?,?,?)"; 
		 try 
		 {
			PreparedStatement ps_new = conn.prepareStatement(sql_new);
			PreparedStatement ps_old = conn.prepareStatement(sql_old);
			for (int j = pTree.size()-1; j >= 0; j--)
			 {
				//System.out.println(j);
				//System.out.println(pNode_Index);
				ps_new.setObject(1, pTree.get(j).get(pNode_Index).st_new);
				ps_new.setInt(2, pTree.get(j).get(pNode_Index).id);
				ps_new.setInt(3, pTree.get(j).get(pNode_Index).fileid);
				ps_new.setString(4, pTree.get(j).get(pNode_Index).time);
				ps_new.setInt(5, pTree.get(j).get(pNode_Index).link_id);
				ps_new.addBatch();
				//ps_new.executeUpdate();
				ps_old.setObject(1, pTree.get(j).get(pNode_Index).st_new);
				ps_old.setInt(2, pTree.get(j).get(pNode_Index).id);
				ps_old.setInt(3, pTree.get(j).get(pNode_Index).fileid);
				ps_old.setString(4, pTree.get(j).get(pNode_Index).time);
				ps_old.setInt(5, pTree.get(j).get(pNode_Index).link_id);
				ps_old.addBatch();
				int beforeid=-1;
				if (j!=0)
				{
					beforeid = pTree.get(j-1).get(0).id;					
				}
				pTree.get(j).get(pNode_Index).beforeid = beforeid;
				tra[pTree.get(j).get(pNode_Index).id] = pTree.get(j).get(pNode_Index);
				pNode_Index = pTree.get(j).get(pNode_Index).root_index;	
			 }
			ps_new.executeBatch();
			ps_old.executeBatch();
			for(int i = 0;i<invalidPoint.size();i++)
			{
				ps_new.setObject(1, invalidPoint.get(i).st_old);
				ps_new.setInt(2, invalidPoint.get(i).id);
				ps_new.setInt(3, invalidPoint.get(i).fileid);
				ps_new.setString(4, invalidPoint.get(i).time);
				ps_new.setInt(5, -1);
				ps_new.addBatch(); 
				ps_old.setObject(1, invalidPoint.get(i).st_old);
				ps_old.setInt(2, invalidPoint.get(i).id);
				ps_old.setInt(3, invalidPoint.get(i).fileid);
				ps_old.setString(4, invalidPoint.get(i).time);
				ps_old.setInt(5, -1);
				ps_old.addBatch(); 
				tra[invalidPoint.get(i).id] = invalidPoint.get(i);
			}
			ps_new.executeBatch();	
			ps_old.executeBatch();
			for (int i = 0; i < adjacentPoint.size(); i++)
			{
				int beforepointid = adjacentPoint.get(i).before_point_id;
				Struct st_link = tra[beforepointid].st_link;
				ClosestPoint cpt = GetCPTAndDist(adjacentPoint.get(i).st_point, st_link); //获取最近点和距离
				ps_new.setObject(1, cpt.struct);
				ps_new.setInt(2, adjacentPoint.get(i).tra_id);
				ps_new.setInt(3, adjacentPoint.get(i).fileid);
				ps_new.setString(4, adjacentPoint.get(i).time);
				ps_new.setInt(5, tra[beforepointid].link_id);
				ps_new.addBatch(); 
				ps_old.setObject(1, adjacentPoint.get(i).st_point);
				ps_old.setInt(2, adjacentPoint.get(i).tra_id);
				ps_old.setInt(3, adjacentPoint.get(i).fileid);
				ps_old.setString(4, adjacentPoint.get(i).time);
				ps_old.setInt(5, tra[beforepointid].link_id);
				ps_old.addBatch(); 
                tra[adjacentPoint.get(i).tra_id]= new Tra_Node(-1, adjacentPoint.get(i).tra_id, tra[beforepointid].link_id, 
                		adjacentPoint.get(i).fileid, adjacentPoint.get(i).time, -1, adjacentPoint.get(i).st_point, 
                		adjacentPoint.get(i).geo_point, cpt.struct, cpt.geo, st_link, tra[beforepointid].geo_link,
                		tra[beforepointid].link_startNodeID, tra[beforepointid].link_endNodeID, tra[beforepointid].oneway, 
                		tra[beforepointid].link_length, -1, null, null, -1);
			}
			ps_new.executeBatch();
			ps_old.executeBatch();
			ps_new.close();
			ps_old.close();
		 } 
		 catch (Exception e)
		 {
		    e.printStackTrace();
		 }
	}	
	
	/**
	 * @param timeExceptionPoint 时间异常点
	 */
	private void SolveTimeExceptionPoint(Hashtable<Integer, List<Tra_Node>> allNodes,List<Tra_Node> timeExceptionPoint,Network wh_roads)
	{		 
		String sql_new = "insert into wh_modifiedTra(new,tra_id,file_id,time,roadid) values(?,?,?,?,?)"; 
		String sql_old = "insert into wh_tras_old(origin,tra_id,file_id,time,road_id) values(?,?,?,?,?)"; 
		 try 
		 {
			PreparedStatement ps_new = conn.prepareStatement(sql_new);
			PreparedStatement ps_old = conn.prepareStatement(sql_old);
	    	SortedSet<Integer> remain = new TreeSet<Integer>();
			for (int i = tra.length-1; i >=0 ; i--) 
			{
				if (tra[i]!=null)
				{
					if (tra[i].beforeid!=-1 && tra[i].beforeid!=i-1)
					{
						//System.out.println("----------------------------------------"+i);
						int beforeid = tra[i].beforeid;
						//*********************************************
						Struct st = null;
						if (tra[i].path!=null)
						{
//							if (tra[i].subPath!=null)
//							{
//								tra[i].subPath.computeGeometry(0.05);
//								JGeometry geo_path = tra[i].subPath.getGeometry();
//								st = JGeometry.storeJS(geo_path,conn);
//							} 
//							else 
//							{
//								tra[i].path.computeGeometry(0.05);
//								JGeometry geo_path = tra[i].path.getGeometry();
//								st = JGeometry.storeJS(geo_path,conn);
//							}
							tra[i].path.computeGeometry(0.05);
							JGeometry geo_path = tra[i].path.getGeometry();
							st = JGeometry.storeJS(geo_path,conn);
						}
						else
						{
							st = tra[i].st_link; //此处需要优化
						}
						//********************************************
						for (int j = beforeid+1; j < i; j++)
						{
							if (tra[j]!=null)
							{
								continue;
							}
							Struct st_old = null;
							JGeometry geo_old = null;
							int fileid=-1;
							String time = null;
							int index = -1;
							for(int m=timeExceptionPoint.size()-1;m>=0;m--)
							{
								if (timeExceptionPoint.get(m).id==j)
								{
									st_old = timeExceptionPoint.get(m).st_old;
									geo_old= timeExceptionPoint.get(m).geo_old;
									fileid = timeExceptionPoint.get(m).fileid;
									time = timeExceptionPoint.get(m).time;
									index = m;
									break;
								}
							}
							if (st_old!=null && index!=-1)
							{
								if (st==null)
								{
									continue;
								}
								ClosestPoint cpt = GetCPTAndDist(st_old, st);							
								timeExceptionPoint.remove(index);
								int linkid = getCloestlinkid(cpt.struct, conn);
								ps_new.setObject(1, cpt.struct);
								ps_new.setInt(2, j);
								ps_new.setInt(3, fileid);
								ps_new.setString(4, time);
								ps_new.setInt(5, linkid);
								ps_new.addBatch();
								ps_new.executeBatch();
								ps_old.setObject(1, st_old);
								ps_old.setInt(2, j);
								ps_old.setInt(3, fileid);
								ps_old.setString(4, time);
								ps_old.setInt(5, linkid);
								ps_old.addBatch();
								ps_old.executeBatch();
								tra[j]=new Tra_Node(-1, j, linkid, fileid, time, -1, st_old, geo_old, cpt.struct, cpt.geo, 
										null, null, -1, -1, null, -1, -1, null, null, -1);
							}
						}
						i=tra[i].beforeid+1;
					}
					else if (i!=0 && tra[i].beforeid ==-1 && tra[i-1]==null )
					{
						int CurrentTraid = i;
						int NextTraid = -1;
						List<Integer> trasid = new ArrayList<Integer>(); //用于记录本次有问题的轨迹点
						for(int n=CurrentTraid-1;n>=0;n--)
						{
							if (tra[n]==null)
							{
								int index = -1;
								for(int m = timeExceptionPoint.size()-1; m >= 0;m--)
								{
									if (timeExceptionPoint.get(m).id==n)
									{
										tra[n]=timeExceptionPoint.get(m);
										index = m;
										break;
									}
								}
								if (index!=-1)
								{
									timeExceptionPoint.remove(index);
									trasid.add(n);
								}
							}
							else if(tra[n]!=null && tra[n].percentage !=-1)
							{
								NextTraid = n;
								break;
							}
						}
						Struct st = null;
						if (NextTraid != -1)
						{
							//*****************************************
							Tra_Node before = tra[NextTraid];
							if (NetworkManager.isReachable(wh_roads, before.link_startNodeID, tra[CurrentTraid].link_endNodeID))
							{
								if (before.link_id == tra[CurrentTraid].link_id)
								{
									st = before.st_link;
								}
								else
								{
									//****************************************************
									if ((before.percentage>0 && before.percentage<1)||(tra[CurrentTraid].percentage>0 && tra[CurrentTraid].percentage<1))
									{
										SubPath subpath = GetSubPath(wh_roads, before, tra[CurrentTraid],5000);
									    if (subpath!=null)
									    {
									    	subpath.computeGeometry(0.05);
										    JGeometry geo = subpath.getGeometry();
											st = JGeometry.storeJS(conn, geo);
										}
									}
									else
									{
										int startNodeID = -1;
										if (before.percentage==0)
										{
											startNodeID = before.link_startNodeID;
										}
										else
										{
											startNodeID = before.link_endNodeID;
										}
										int endNodeID = -1;
										if (tra[CurrentTraid].percentage == 0)
										{
											endNodeID = tra[CurrentTraid].link_startNodeID;
										}
										else 
										{
											endNodeID = tra[CurrentTraid].link_endNodeID;
										}
										if (startNodeID != -1 && endNodeID != -1) 
										{
											if (startNodeID != endNodeID)
											{
												Path path = NetworkManager.shortestPath(wh_roads, startNodeID, endNodeID);
											    if (path!=null)
											    {
											    	path.computeGeometry(0.05);
												    JGeometry geo = path.getGeometry();
													st = JGeometry.storeJS(conn, geo);
												}
											}										
										}
									}
									
								}
							}
							if (st == null)
							{
								for(int q=NextTraid+1;q<CurrentTraid;q++)
								{
									remain.add(q);
								}
								i=NextTraid+1;
								continue;
							}
							for(int p=NextTraid+1;p<CurrentTraid;p++)
							{
								if (tra[p].link_id!=-1 || !trasid.contains(p))
								{
									continue;
								}
								ClosestPoint cpt = GetCPTAndDist(tra[p].st_old, st);	
								int linkid = getCloestlinkid(cpt.struct, conn);
								ps_new.setObject(1, cpt.struct);
								ps_new.setInt(2, tra[p].id);
								ps_new.setInt(3, tra[p].fileid);
								ps_new.setString(4, tra[p].time);
								ps_new.setInt(5, linkid);
								ps_new.addBatch();
								ps_old.setObject(1, tra[p].st_old);
								ps_old.setInt(2, tra[p].id);
								ps_old.setInt(3, tra[p].fileid);
								ps_old.setString(4, tra[p].time);
								ps_old.setInt(5, linkid);
								ps_old.addBatch();
							}
							ps_new.executeBatch();
							ps_old.executeBatch();
						}
						i=NextTraid+1;
					}
				}
				else
				{
					//*********************************
					remain.add(i);
				}								
			}
			//**************************************
			Object[] remainpoint = remain.toArray();
			for(int n= 0;n<remain.size();n++)
			{
				int key = (int)(remainpoint[n]);
				List<Tra_Node> gpsNodes = new ArrayList<Tra_Node>();
				gpsNodes = allNodes.get(key);
			    //获取前一结点
			    Tra_Node before = null;
			    for (int j = key-1; j >=0; j--)
			    {
					if (tra[j].link_id!=-1)
					{
						before = tra[j];
						break;
					}
				}
			    double cost = 0;
			    Tra_Node suitableNode = null;
			    if (before!=null)
			    {
					for (int j = 0; j < gpsNodes.size(); j++) 
					{
						SubPath subpath = null; //最短路径
						double dist_route = 0;
						//double dist_Adjacent = before.geo_old.distance(gpsNodes.get(j).geo_old, 0.05, "true");
						double dist_Adjacent = getDistance(before.st_old, gpsNodes.get(j).st_old, conn);
						if (NetworkManager.isReachable(wh_roads, before.link_startNodeID, gpsNodes.get(j).link_endNodeID))
						{
							if (before.link_id == gpsNodes.get(j).link_id)
							{
								dist_route = Math.abs(gpsNodes.get(j).percentage-before.percentage) * before.link_length;
							}
							else
							{
								if ((before.percentage>0 && before.percentage<1)||(gpsNodes.get(j).percentage>0 && gpsNodes.get(j).percentage<1))
								{
									subpath = GetSubPath(wh_roads, before, gpsNodes.get(j),dist_Adjacent*5);
								}
								else
								{
									int startNodeID = -1;
									if (before.percentage==0)
									{
										startNodeID = before.link_startNodeID;
									}
									else
									{
										startNodeID = before.link_endNodeID;
									}
									int endNodeID = -1;
									if (gpsNodes.get(j).percentage == 0)
									{
										endNodeID = gpsNodes.get(j).link_startNodeID;
									}
									else 
									{
										endNodeID = gpsNodes.get(j).link_endNodeID;
									}
									if (startNodeID != -1 && endNodeID != -1) 
									{
										if (startNodeID == endNodeID)
										{
											dist_route = 0;
										}
										else 
										{
											Path path = NetworkManager.shortestPath(wh_roads, startNodeID, endNodeID);
											if (path!=null)
											{
												dist_route = path.getCost();
											}
										}										
									}
								}
							}
						}
						//求路径距离
						if (subpath != null)
						{
							dist_route = subpath.getCost();
						}
						//路径距离为0 
						if (dist_route == 0)
						{
							continue;
						}			
						double dt = Math.abs(dist_Adjacent-dist_route);
						double tp_cost = 0.039 * Math.exp(-dt * 0.039) * constvalue; //转移概率
						if (gpsNodes.get(j).cost * tp_cost* before.cost >cost) 
						{
							cost = gpsNodes.get(j).cost * tp_cost* before.cost;
							if (subpath != null)
							{
								gpsNodes.get(j).path = subpath.getReferencePath();
							}
							else
							{
								gpsNodes.get(j).path = null;
							}							
							gpsNodes.get(j).subPath = subpath;
							suitableNode = gpsNodes.get(j);
						}			
					}
				}
			    if (suitableNode == null)
			    {
			    	if (gpsNodes.size() != 0)
			    	{
			    		ps_new.setObject(1, gpsNodes.get(0).st_old);
						ps_new.setInt(2, key);
						ps_new.setInt(3,  gpsNodes.get(0).fileid);
						ps_new.setString(4, gpsNodes.get(0).time);
						ps_new.setInt(5, -1);
						ps_new.addBatch();
						ps_new.executeBatch();
						ps_old.setObject(1, gpsNodes.get(0).st_old);
						ps_old.setInt(2, key);
						ps_old.setInt(3, gpsNodes.get(0).fileid);
						ps_old.setString(4, gpsNodes.get(0).time);
						ps_old.setInt(5, -1);
						ps_old.addBatch();
						ps_old.executeBatch();
					}
				}
			    else
			    {
			    	tra[key] = suitableNode;
				    ps_new.setObject(1, suitableNode.st_new);
					ps_new.setInt(2, key);
					ps_new.setInt(3, suitableNode.fileid);
					ps_new.setString(4, suitableNode.time);
					ps_new.setInt(5, suitableNode.link_id);
					ps_new.addBatch();
					ps_new.executeBatch();
					ps_old.setObject(1, suitableNode.st_old);
					ps_old.setInt(2, key);
					ps_old.setInt(3, suitableNode.fileid);
					ps_old.setString(4, suitableNode.time);
					ps_old.setInt(5, suitableNode.link_id);
					ps_old.addBatch();
					ps_old.executeBatch();
			    }
			   
			}
			ps_new.close();
			ps_old.close();
		 } 
		 catch (Exception e)
		 {
		    e.printStackTrace();
		    logger.error(Arrays.toString(e.getStackTrace()));
		 }
	}
	
	   private int getCloestlinkid(Struct tra,Connection conn)
	   {	   
			try
			{
				CallableStatement cs = conn.prepareCall( "{ call getclosestlinkid2(?,?) }" );
				cs.setObject(1, tra);
				cs.registerOutParameter ( 2, OracleTypes.NUMBER);
				cs.execute ();
				int linkid = cs.getInt(2);
				cs.clearBatch();
				cs.close();
				return linkid;
			} 
			catch (SQLException e) 
			{
				//System.out.println("求最近link id出现异常："+e.getMessage());
				return -1;
			}
	   }
	
	//距离相差太近点
	private class AdjacnetPoint
	{
		public int before_point_id;
		public Struct st_point;
		public JGeometry geo_point;
		public int tra_id;
		public int fileid;
		public String time;
		/**
		 * @param _before_point_id 距离最近的有效点，一般为前一点
		 * @param _st_point 本点的Struct
		 * @param _geo_point 本点的geom
		 * @param _tra_id 本点id
		 * @param _fileid 本点文件id
		 * @param _time 本点采集时间
		 */
		public AdjacnetPoint(int _before_point_id,Struct _st_point,JGeometry _geo_point,int _tra_id,int _fileid,String _time)
		{
			before_point_id = _before_point_id;
			st_point = _st_point;
			tra_id = _tra_id;
			fileid = _fileid;
			time = _time;
			geo_point = _geo_point;
		}
	}
	/**
	 * @param networkname 网络数据集名称
	 * @param link_id  link编号
	 * @param pt link上的一点
	 * @return 从起点到pt所占的百分比
	 */
	private double GetPercentage(String networkname,int link_id,Struct pt)
	{
		double percentage = 0;
		try
		{
			CallableStatement cs = conn.prepareCall( "{ call getpercentage(?,?,?,?) }" );
			cs.setString(1, networkname);
			cs.setObject(2, pt);
			cs.setInt(3, link_id);
			cs.registerOutParameter ( 4 , OracleTypes.NUMBER);
			cs.execute ();
			percentage = cs.getDouble(4);
			cs.clearBatch();
			cs.close();
			return percentage;
		} 
		catch (SQLException e) 
		{
			System.out.println("求百分比出现异常："+"link_id为" + link_id +e.getMessage());
			logger.debug("求百分比出现异常："+"link_id为" + link_id +e.getMessage());
			return percentage;
		}
	}
	
	/**
	 * @param time 时间字段
	 * @return 返回标准格式时间
	 */
	private Date GetDate(String time)
	{
		Date date = new Date();
        SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd hh:mm:ss");
        try 
        {
			date = sdf.parse(time);
			return date;
		}
        catch (ParseException e) 
        {
        	System.out.println("时间转换失败！");
        	logger.debug("时间转换失败！");
			return null;
		}
	}
	
	/**
	 * @param struct_point 轨迹点
	 * @param struct_line  路段
	 * @return 距离和最近点
	 */
	private ClosestPoint GetCPTAndDist(Struct struct_point,Struct struct_line)
	{
		ClosestPointAndDist CPTAndDist = new ClosestPointAndDist(conn, struct_point, struct_line);
	    return CPTAndDist.CPTAndDist();
	}		
	
	private double getDistance(Struct p1,Struct p2,Connection conn)
	   {
		   double dis = 0;
			try
			{
				CallableStatement cs = conn.prepareCall( "{ call getdistance(?,?,?) }" );
				cs.setObject(1, p1);
				cs.setObject(2, p2);
				cs.registerOutParameter ( 3, OracleTypes.NUMBER);
				cs.execute ();
				dis = cs.getDouble(3);
				cs.clearBatch();
				cs.close();
				return dis;
			} 
			catch (SQLException e) 
			{
				System.out.println("求距离出现异常："+e.getMessage());
				logger.debug("求距离出现异常："+e.getMessage());
				return dis;
			}
	   }
}
