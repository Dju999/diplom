#-*-coding:utf8-*-
import _mysql#загруженный пакет MySQLdb
import csv#чтение csv
#набор пакетов для кластеризации и построения графики
from pylab import plot,show
from numpy import vstack,array
from numpy.random import rand
from scipy.cluster.vq import kmeans,vq
#для гистограммы
import matplotlib.pyplot as plt
from matplotlib import rc#для кирилицы в подписях
from matplotlib.pyplot import legend
#работа с файловой системой
import os
#для логарифмирования
import math
#генерация нормальных СВ
from scipy.stats import norm
from scipy.stats import normaltest
from scipy.stats import skew
from scipy.stats import kurtosis
from scipy.stats import chisquare
from scipy.stats import chi2
import random

rc('font',**{'family':'serif', 'size'   : 18})
rc('text', usetex=True)
rc('text.latex',unicode=True)
rc('text.latex',preamble='\usepackage[utf8]{inputenc}')
rc('text.latex',preamble='\usepackage[russian]{babel}')

def pearson_test(sample,groups, mu, sigma,alpha):
	sample=list(map(lambda x: round(x,2),sample))
	#print groups
	result = list()
	n=len(sample)
	print "Объём выборки %d"%n
	k=(len(groups)-1)-1#число степеней свободы
	sample.sort()
	#print sample
	count = 0
	print groups
	for i in groups:
		p_i=norm.cdf(sample[count+i-1],0,sigma)-norm.cdf(sample[count],0,sigma)
		#print "Интервал:[%.2f;%.2f],%d элементов"%(sample[count],sample[count+i-1],i)
		result.append((i-n*p_i)**2/(n*p_i))
		count=count+i
	P = sum(result)
	q = chi2.ppf(alpha,k)
	print "Значение статитики %.2f, значение квантили %.2f на уровне значимости a=%.2f => гипотеза о нормальности %s."%(P,q,alpha,("принимается" if P<q else "не принимается"))
	return P

class IRT:
	def __init__(self, dbparams):
		self.dbprm = dbparams
		self.item_list = list()
		self.pseudo_users = list()
		self.users = list()
		self.ind = list()
	def get_db_params(self):
		return self.dbprm
	def dbconnection(self):
		self.db=_mysql.connect(host=self.dbprm["host"],user=self.dbprm["user"], passwd=self.dbprm["passwd"],db=self.dbprm["db"])
	def fit_item(self,prefix):#отбираем задачи, время ответа на которые распределено нормально
		f=open('table','w')
		#print round(chi2.ppf(0.95,2),3)#значение квантили
		fit = list()
		self.get_item_list(prefix)#получаем список доступных задач
		task_list = self.item_list
		for item in task_list:
			test_data = list(map(lambda x: math.log(x),self.get_irt_for_item(prefix,item)[0]))
			if len(test_data)>20 and len(test_data)<100:
				try:
					s = round(skew(test_data),2)
					k = round(kurtosis(test_data,fisher=False),2)
					Z = normaltest(test_data)
					C = chisquare(test_data)
					if Z[1]>0.77:
						#print "Item: %s;\tStat: %s;\tp-value: %s;\tchi-square: %s;\tp-value: %s;\tskew: %s;\tkurtosis: %s"%(item,round(Z[0],3),round(Z[1],2),round(C[0],2),round(C[1],2),s,k)
						#f.write("\hline\n%s & %s & %s \\\\ \n"%(item,round(Z[0],3),round(Z[1],2)))
						f.write("\hline\n%s & %s & %s & %s & %s \\\\ \n"%(item,len(test_data),s,k,round(Z[0],3)))
						fit.append(item)
				except ValueError:
					pass
		f.write("\hline");f.close()
		return fit
	def fit_validation(self,task_list_origin,irt_table, beta, tau, sigma, mu):#из списка задач выкидываем ту, которая сильнее всего "портит" гипотезу о нормальности
		evil1=666#значение статистики
		task_list=task_list_origin[:]
		for i in xrange(len(task_list_origin)):
			print "Задача %d,%s"%(i,task_list[i])
			del task_list[i]
			irt_table = self.get_irt_table('matan',task_list)
			mu = self.hat_mu(irt_table)
			beta = self.hat_beta(irt_table,mu)
			tau = self.hat_tau(irt_table,mu)
			sigma = self.hat_sigma(irt_table, tau, beta, mu)
			evil2 = self.check_fit(irt_table, beta, tau, task_list, sigma, mu)
			evil1 = (evil1 if evil1<evil2 else evil2)
			task_list=task_list_origin[:]
		return evil1
	def get_irt_table(self, prefix, item_list):
		data = {}.fromkeys(item_list,0)
		self.dbconnection()
		for item in item_list:
			tmp={}
			#если юзер ответил неск. раз - усредняем значение
			self.db.query("SELECT AVG(irt), user FROM "+prefix+"_www_irt WHERE item='"+item+"' GROUP BY user ORDER BY user")
			r=self.db.use_result()
			line = r.fetch_row()
			while line:
				if float(line[0][0]) < 3600.0:
					tmp.update({int(line[0][1]): math.log(float(line[0][0]))})
				line = r.fetch_row()
			data[item] = tmp
		return data
	def hat_mu(self, data):
		res = list()
		for key1, i in data.items():
			for key2, j in i.items():
				res.append(j)
		return sum(res)/len(res)
	def hat_beta(self, data, mu_hat):
		keys = data.keys()
		res = {}.fromkeys(keys,0)#совпадает с количеством задач
		for key,item in data.items():
			tmp = item.values()
			res[key] = round(sum(tmp)/len(tmp)-mu_hat,2)
		return res
	def hat_tau(self, data, mu_hat):
		res = {}#совпадает с количеством студентов
		for key,item in data.items():
			tmp = item.keys()
			for user in tmp:
				if not res.has_key(user):
					res.update({user:list()})
		for key1, item in data.items():
			for user, ir_time in item.items():
				res[user].append(ir_time)
		for key, ir_time in res.items():
			res.update({key: sum(ir_time)/len(ir_time) - mu_hat})
		return res
	def hat_sigma(self, data, tau_hat, beta_hat, mu_hat=0):
		res = list()
		for item, i in data.items():
			for user, j in i.items():
				res.append((j - tau_hat[user] - beta_hat[item]-mu_hat)**2)
		return math.sqrt(sum(res)/len(res))
	def check_fit(self, data, beta, tau, task_list, sigma, hat_mu, alpha):
		f = open('errors.txt','w')
		errors = list()
		for item, i in data.items():
			for user, j in i.items():
				errors.append(j - tau[user] - beta[item]-hat_mu)
				f.write("%f "%errors[-1])
		#Z = normaltest(errors)
		groups=self.irt_gist(errors,"error.png","значение $\\varepsilon$","частота")#получаем кол-во элементов в каждом столбце
		#print Z
		f.close()
		P = pearson_test(errors,groups, hat_mu, sigma, alpha)
		return P
	def latex_table(self, hat_mu, beta, tau, task_list, sigma, L, stud_id, forecast, model):
		f = open('latex_table.txt','w')
		f.write("id студента: %s\s"%stud_id)
		f.write("временной параметр студента: %s\s"%tau[stud_id])
		f.write("параметр компроментации: %s\n"%L)
		f.write("Задача №&    Распределение времени ответа &            Прогнозное время &           Фактическое время                       & Доверит. интервал  &                         Вывод\\\\ \n")
		f.write("\\hline\n")
		for ind, i in enumerate(task_list):
			f.write("%d &           N(%.2f, %.2f) &                             %.2f &                         %.2f  &                            $[%.2f;+\infty]$ &                                 %s\\\\ \n"%
				   (ind + 1,(hat_mu+beta[i]+tau[stud_id]),sigma, forecast[task_list.index(i)], model[task_list.index(i)]*(-1), norm.ppf(0.05)*sigma+hat_mu+beta[i]+tau[stud_id],  ("+" if model[task_list.index(i)]*(-1) > norm.ppf(0.05)*sigma+hat_mu+beta[i]+tau[stud_id] else "-")     ))
			f.write("\\hline\n")
		f.close()
	def print_table(self, data, type_of_table, tau=None, beta=None, number_of_stud=10, number_of_items=10):
		'''
			рандомно выбираем заданное количество студентов/зачач
		'''
		rdict = {}
		f = open(type_of_table+'.txt','w')
		f.write('{|c')
		if type_of_table == "student_params":
			for i in xrange(number_of_stud):#записываем в файл первую строчку и одновременно генерим рандомный список студентов
				sex = True; ind = None;
				while sex:
					ind = random.choice(tau.keys())
					if ind not in rdict:
						rdict.update({ind:tau[ind]})
						sex = False
				f.write('\t|c')
			f.write('|}\n\\hline\n')
			#заполняем первую строку
			f.write('id студента  \t')
			for key, i in rdict.items():
				f.write("\t&\t%s\t"%key)
			f.write('\\\\ \n\\hline\n$\\tau$ студента\t')
			#заполняем вторую строку
			for key, i in rdict.items():
				f.write("\t&\t%s\t"%round(i,2))
		if type_of_table == "item_params":
			for i in xrange(number_of_items):#записываем в файл первую строчку и одновременно генерим рандомный список студентов
				sex = True; ind = None;
				while sex:
					ind = random.choice(beta.keys())
					if ind not in rdict:
						rdict.update({ind:beta[ind]})
						sex = False
				f.write('\t|c')
			f.write('|}\n\\hline\n')
			#заполняем первую строку
			f.write('Задача  \t')
			for key, i in rdict.items():
				f.write("\t&\t%s\t"%key)
			f.write('\\\\ \n\\hline\n$\\beta$ задачи\t')
			#заполняем вторую строку
			for key, i in rdict.items():
				f.write("\t&\t%s\t"%round(i,2))
		if type_of_table == "irt_params":
			for stud, i in tau:
				for item, j in beta:
					pass
		f.write("\\\\\n\\hline")
		f.close()
		return rdict
	def irt_gist(self,data,name,x_label_name,y_label_name):
		try:
			plt.hist(data,bins=8, color='#BCF5BC', normed=True)
			#print n.tolist()
			#print bins.tolist()
			plt.xlabel(x_label_name)
			#plt.ylabel(y_label_name)
			plt.savefig(name, dpi = 50)
			plt.clf()
			plt.cla()
			n, bins, patches = plt.hist(data,bins=8)#получаем разбиение
			return n.tolist()
		except IndexError:
			pass
	def plot_forecast(self,task_list,forecasted,modeled, diff):
		ab = [x*1.6 for x in xrange(1,len(task_list)+1)]#координатная ось, расширение с коэфф. 1.3 - чтобы всё влезло
		#прогноз - в верхней полуплоскости, наблюдаемое значение в нижней(в таком формате возвращается из real_test)
		plt.bar(ab, forecasted, facecolor='#9999ff', edgecolor='white', align = 'center')
		plt.bar(ab, modeled,    facecolor='#ff9999', edgecolor='white', align = 'center')
		for x, y in zip(ab, forecasted):
			plt.text(x+0.0, y+1.5, '%.2f' % y, ha='center', va='top', fontsize=12)#верхняя часть
		for x, y in zip(ab, modeled):
			plt.text(x+0.0, y-1.0, '%.2f' % y, ha='center', va='top', fontsize=12)#нижняя часть
		plt.xticks(ab,task_list)
		plt.xlabel('Задачи')
		plt.ylabel('$\ln t_{ij}$')
		plt.ylim(-23.25, +23.25)
		plt.savefig("forecast.eps", dpi = 50)
		plt.clf()
		plt.cla()
		#теперь отобразим в одной полуплоскости
		modeled = [(-1)*x for x in modeled]
		a1 = [x-0.2 for x in ab]#столбики теперь размещаются не в одной точке, а рядом друг с другом
		a2 = [0.2+x for x in ab]#столбики теперь размещаются не в одной точке, а рядом друг с другом
		plt.bar(a1, forecasted, facecolor='#9999ff', edgecolor='white', align = 'center',width=0.4, label="$t_{ij}^{\mbox{п}}$ прогноз")
		plt.bar(a2, modeled,    facecolor='#ff9999', edgecolor='white', align = 'center',width=0.4, label="$t_{ij}^{\mbox{м}}$ наблюдения")
		#размещаем подписи - над каждым столбцом
		for x, y in zip(a1, forecasted):
			plt.text(x-0.2, y+0.5, '%.2f' % y, ha='center', va='top', fontsize=12)
		for x, y in zip(a2, modeled):
			plt.text(x+0.2, y+0.5, '%.2f' % y, ha='center', va='top', fontsize=12)
		plt.xticks(ab,task_list)
		plt.xlabel('Задачи')
		plt.ylabel('$\ln t_{ij}$')
		plt.ylim(-0.25, +12.00)
		legend(bbox_to_anchor=(0.6, 0.85), loc=3, borderaxespad=0.,prop={'size':16})
		plt.savefig("forecast-2.eps", dpi = 50)
		plt.clf()
		plt.cla()
		#разница между моделированной и предсказанной величиной
		plt.bar(ab, diff, facecolor='#2EDC91', edgecolor='white', align = 'center')
		for x, y in zip(ab, diff):
			plt.text(x+0.0, y+0.25, '%.2f' % y, ha='center', va='top', fontsize=12)
		plt.xticks(ab,task_list)
		plt.xlabel('Задачи')
		plt.ylabel('Ошибка прогноза')
		plt.ylim(-6.25, +6.25)
		plt.savefig("fdiffer.eps", dpi = 50)
		plt.clf()
		plt.cla()
	def save_gist(self, prefix):
		save_dir = "./GIST/"+prefix+"/"#директория для сохранения
		self.get_item_list(prefix)#получаем список задач по предмету prefix
		for i in self.item_list:#проходимся по всем задачам, сохраняем картинки
			data = self.get_irt_for_item(prefix,i)
			if len(data[0])>15:
				#переходим к логарифмической шкале времени
				try:
					log_data = list(map(lambda x: math.log(x),data[0]))
				except ValueError:#если попалось нулевое значение - просто выкидываем из списка
					data = list(filter(lambda x: x!=0,data[0]))
					log_data = list(map(lambda x: math.log(x),data))
				self.irt_gist(log_data, save_dir+i+".png")
	def get_irt_for_item(self, prefix, item):#irt для одной конкретной задачи
		r1 = list()# все ответы
		r2 = list()#только неправильные
		self.dbconnection()
		self.db.query("SELECT irt FROM "+prefix+"_www_irt WHERE item='"+item+"' AND answer=0 ORDER BY irt")#вытаскиваем только неправильные
		r=self.db.use_result()
		line = r.fetch_row()
		while line:
			if int(line[0][0]) < 1300000000:
				r1.append(int(line[0][0]))
			line = r.fetch_row()
		self.db.query("SELECT irt FROM "+prefix+"_www_irt WHERE item='"+item+"' ORDER BY irt")#вытаскиваем все вместе
		r=self.db.use_result()
		line = r.fetch_row()
		while line:
			if int(line[0][0]) < 1300000000:
				r2.append(int(line[0][0]))
			line = r.fetch_row()
		return (r2, r1)#r1 - неправильные, r2 - все вместе
	def get_item_list(self,prefix):
		self.item_list = list()
		self.dbconnection()
		self.db.query("SELECT DISTINCT item FROM "+prefix+"_www_irt ORDER BY item");
		r=self.db.use_result()
		line = r.fetch_row()
		while line:
			self.item_list.append(line[0][0])
			line = r.fetch_row()
	'''
	def get_stat(self, prefix):#получаем статистику правильных и неправильных ответов по всем задачам предмета
		self.get_item_list(prefix)
		for i in self.item_list:
			self.get_irt_for_item(prefix, i)
	def get_irt_by_item(self, prefix):#irt для всех доступных задач
		self.get_item_list(prefix)#получаем список доступных задач
		for i in self.item_list:#проходимся по списку задач
			self.dbconnection()
			self.db.query("SELECT irt FROM "+prefix+"_www_irt WHERE item='"+i+"' ORDER BY irt")#время, которое разные студенты затратили на задачу, выбрасываем корректные ответы
			r=self.db.use_result()
			line = r.fetch_row()
			f=open('./ITEMS/'+prefix+'/'+i+'.txt', 'w')
			while line:#записываем в соответствующий файл время ответа на задачу (в порядке возрастания)
				if int(line[0][0]) < 1300000000:
					f.write(line[0][0]+'\n')
				line = r.fetch_row()
			f.close()
	def get_raw_data(self, prefix):#импорт из БД в файл .csv
		#открываем файл для записи
		result_file = open('./CSV/'+prefix+'.csv', 'w')
		result_file.write('rec_id,item_id,stud_id,result,irt\n')
		self.dbconnection()
		self.db.query("SELECT * FROM "+prefix+"_www_irt ORDER BY item");
		r=self.db.use_result()
		line = r.fetch_row()
		countall = 0
		count = 0
		while line:
			if int(line[0][4]) > 1300000000:
				count = count+1
			else:
				for i in line[0]:
					result_file.write(i)
					if i!=line[0][4]:
						result_file.write(",")
					else:
						result_file.write("\n")
			countall = countall+1
			line = r.fetch_row()
		result_file.close()
	def get_pseudo_users(self, ts, step):
		self.dbconnection()
		#запрос к sql:  select distinct user from matan_www_irt where item='4.2.1' or item='5.2.6' ...;
		qu = "select distinct user from matan_www_irt where"
		for i in xrange(len(ts)):
			if i == 0:
				qu = qu+" item='"
			else:
				qu = qu+" or item='"
			qu = qu+ts[i]+"'"
		qu = qu+";"
		self.db.query(qu)#вытаскиваем пользователей
		r=self.db.use_result()
		line = r.fetch_row()
		while line:
			self.users.append(line[0][0])
			line = r.fetch_row()
		for i in xrange(0,len(self.users),step):
			index = int(round(rand()*1000,0))#рандомно генерим индекс псевдопользователя
			self.ind.append(index)
		return 0
	def test(self):#тестируем всякое
		task_list = ['4.2.1', '5.2.6', '6.2.5', '6.2.8', '6.2.10', '8.1.5', '8.2.1', '8.2.3', '8.2.5', '8.2.6', '9.2.2', '9.2.1', '9.2.3', '9.2.9', '10.2.1', '10.2.2']
		c = list()
		a = {'11':[0,1,2],'12':[2,3,4],'13': [5,6,7,8,9]}
		b = {'11':[10,11,12,13,14],'12': [15,16,17,18],'13': [19,20,21]}
		c.append(a)
		c.append(b)
		task_list = ['4.2.1', '5.2.6']
		time_list1 = {}.fromkeys(task_list,0)#создаём словарь с заданными ключами и заполняем значениями по дефолту
		time_list2 = {}.fromkeys(a.keys(),0)#создаём словарь с заданными ключами и заполняем значениями по дефолту
		M = len(c[0])#количество студентов
		k = len(c)#количество задач
		M = len(c[0])#количество студентов
		for i in c:
			tmp = list()
			for key, j in i.items():
				for k in j:
					tmp.append(k)
			print tmp
			time_list1[task_list[c.index(i)]]=sum(tmp)/M
	def get_irt(self, ts, step, item):#
		self.dbconnection()
		#получаем список "псевдостудентов" и время ответа для каждого из них
		res = self.users
		result = {}.fromkeys(self.ind)#словарь, в котором храним результат
		for i in xrange(0,len(res),step):
			tmp = list()#сохраняем вемя ответа
			if i<(len(res)-step):
				for j in xrange(step):
					qq = "SELECT irt FROM matan_www_irt WHERE user=%s AND item='%s'" % (res[i+j], item)
					self.db.query(qq)
					rr=self.db.use_result()
					line = rr.fetch_row()
					while line:
						if int(line[0][0]) < 1300000000:
							tmp.append(math.log(float(line[0][0])))
						line = rr.fetch_row()
			else:
				for k in xrange(i,len(res)):
					qq = "SELECT irt FROM matan_www_irt WHERE user=%s AND item='%s'" % (res[k], item)
					self.db.query(qq)
					rr=self.db.use_result()
					line = rr.fetch_row()
					while line:
						if int(line[0][0]) < 1300000000:
							tmp.append(math.log(float(line[0][0])))
						line = rr.fetch_row()
			result[self.ind[i/step]] = tmp
		return result	
	'''
	def clean_dir(self,dirname,subdir):
		for i in subdir:
			filedir = "./"+dirname+"/"+i
			files = os.listdir(filedir)
			for j in files:
				os.remove(filedir+"/"+j)
	def forecast(self, hat_mu, beta, tau, sigma,L=0):
		#print hat_mu+beta+tau+L
		return round(norm.rvs(loc=(hat_mu+beta+tau+L), scale=sigma, size=1)[0],2)
	def real_test(self, hat_mu, beta, tau, task_list, sigma, L, stud_id):
		modeled = list()
		forecasted = list()
		for i in task_list:
			if task_list.index(i)==2 or task_list.index(i)==5 or task_list.index(i)==8:
				modeled.append(self.forecast(hat_mu, beta[i], tau[stud_id], sigma, L)*(-1))
			else:
				modeled.append(self.forecast(hat_mu, beta[i], tau[stud_id], sigma)*(-1))
			forecasted.append(self.forecast(hat_mu, beta[i], tau[stud_id], sigma))
		differ = list(map(lambda x,y: round((-1)*x-y,2), modeled,forecasted))
		return (modeled,forecasted,differ)
	def save_params(self, hat_mu, beta, tau, task_list, sigma, L, forecast, model):
		f=open('result.txt', 'w')
		f.write("pseudo_users\t=\t[")
		for i in self.ind:
			f.write("%s,\t"%i)
		f.write("]\n")		
		f.write("\mu\t=\t%s\n"%hat_mu)
		f.write("\sigma\t=\t%s\n"%sigma)
		f.write("L\t=\t%s\n"%L)
		f.write("beta\t=\t[")
		for i in beta:
			f.write("%s,\t"%i)
		f.write("]\n")
		f.write("tau\t=\t[")
		for i in tau:
			f.write("%s,\t"%i)
		f.write("]\n")
		f.write("forecast\t=\t[")
		for i in forecast:
			f.write("%s,\t"%i)
		f.write("model\t=\t[")
		f.write("]\n")
		for i in model:
			f.write("%s,\t"%i)
		f.write("]\n")
		for i in list(map(lambda x,y: x-y, forecast, model)):
			f.write("%s,\t"%i)
		f.write("]\n")
		f.close


dbparams = {"host": "localhost","user": "pythonuser","passwd": "pythonpass","db": "distance"}
irt = IRT(dbparams)#инициализация класса
task_list = irt.fit_item('matan');del task_list[4]#выбираем задачи с высоким p-value
irt_table = irt.get_irt_table('matan',task_list)
mu = irt.hat_mu(irt_table)
beta = irt.hat_beta(irt_table,mu)
tau = irt.hat_tau(irt_table,mu)
sigma = irt.hat_sigma(irt_table, tau, beta, mu)
irt.check_fit(irt_table, beta, tau, task_list, sigma, mu, 0.95)



subj = ["twims","linal","matan"]#список предметов, по которым хотим выгрузить данные из MySQL в csv
#for i in subj:
	#irt.get_raw_data(i)
	#irt.get_irt_by_item(i)
	#irt.save_gist(i)
	#pass
#irt.clean_dir("ITEMS",subj)
#irt.get_stat("matan")
#arr = irt.get_irt_for_item('matan','4.2.9')
#print arr[0]
#data = list(map(lambda x: math.log(x),arr[0]))
#print arr[0]
#irt.irt_gist(data,"4.2.9.png")
#irt.get_vis('matan','6.2.3',2)
#irt.save_gist("matan")

"""
alpha = 0.9#нормировочный коэффициент
sigma = alpha*irt.hat_sigma(irt_table,mu, tau, beta)/len(task_list)
print mu
print beta
print tau
"""
"""
#для случаев множественного прогноза(для неск. студентов)
if False:
	student_params = irt.print_table(irt_table, "student_params", tau, number_of_stud=20)
item_params = irt.print_table(irt_table, "item_params", beta=beta, number_of_items=9)
#результаты эксперимента
stud_id = random.choice(tau.keys())#выбираем случайного студента
L = (-1)*sigma*6.0
task_list = item_params.keys()#формируем список задач, которые попали в тест
#данные для диаграмм
r = irt.real_test(mu, beta, tau, task_list, sigma, L, stud_id)
model = r[0]
forecast = r[1]
differ = r[2]
#сохраняем диаграммы
irt.plot_forecast(task_list,forecast,model, differ)
irt.save_params(mu, beta, tau, task_list, sigma, L, forecast, model)
irt.latex_table(mu, beta, tau, task_list, sigma, L, stud_id, forecast, model)
"""
""" 
sigma = irt.hat_sigma(irt_table, mu, tau, beta)
print sigma
"""