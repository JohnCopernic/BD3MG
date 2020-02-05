
p:']�5  �               @   s�   d  d l  Z d  d l m Z m Z m Z m Z d  d l m Z d  d l Z d  d l m	 Z	 d  d l
 m
 Z
 d  d l Z e j j d � d Z e
 e j d d d	 g � e j d d d
 g � � Z Gd d �  d � Z d S)�    N)�Voperator3D�Vhoperatoradj3D�Vtoperatoradj3D�Vvoperatoradj3D)�	Critere3D)�convolve)�
Gaussian3DzB/home/mathieuchalvidal/PycharmProjects/untitled/BP3MG/FlyBrain.mat�I�   �
   �   c               @   s@   e  Z d  Z d d �  Z d d �  Z d d �  Z d d �  Z d	 S)
� MajorizeMinimizeMemoryGradient3Dc             C   s�  | |  _  | |  _ | |  _ | |  _ | |  _ | |  _ | |  _ | |  _ |	 |  _ |
 |  _	 | |  _
 | |  _ | |  _ | |  _ | |  _ g  |  _ g  |  _ g  |  _ g  |  _ g  |  _ g  |  _ g  |  _ |	 j \ |  _ |  _ |  _ d |  _ d |  _ t d � t d � t d � | d } | d } | d k r>t d � n� | d k rWt d	 � nz | d
 k rpt d � na | d k r�t d � nH | d k r�t d � n/ | d k r�t d � n | d k r�t d � | d k r�t d � na | d k rt d � nH | d
 k rt d � n/ | d k r5t d � n | d k rKt d � t d | � t d | � t d | � t d | � t d | d | � d S) a�  
        :param y: observed data
        :param H: Gaussian blur operator
        :param H_adj: adjoint of the gaussian blur operator
        :param eta: regularization parameter on the distance to the hypercube [xmin, xmax]
        :param lambda_: regularization parameter on the horizontal total variation norm
        :param delta: second regularization parameter on the horizontal total variation norm
        :param kappa: regularization parameter on the vertical total variation norm
        :param phi: choice of regularization function
        :param x: estimation
        :param xstar: ground truth
        :param xbar: minimal estimation
        :param xmin: lower bound on image pixel value
        :param xmax: upper bound on image pixel value
        :param NbIt: nb of iteration for the algorithm to reach
        :param timemax: maximal time of computation
        :return:
        g:�0�yE>�   z(****************************************z+Majorize-Minimize Memory Gradient Algorithmz-> STANDARD VERSION <-r   z%phixy(u) =  (1-exp(-u^2/(2*delta^2)))�   z"phixy(u) = (u^2)/(2*delta^2 + u^2)�   z#phixy(u) = log(1 + (u^2)/(delta^2))�   z#phixy(u) =  sqrt(1 + u^2/delta^2)-1r
   zphixy(u) = 1/2 u^2�   z2phixy(u) = 1-exp(-((1 + (u^2)/(delta^2))^(1/2)-1))�   z%phixy(u) =  (1 + u^2/delta^2)^(1/4)-1z$phiz(u) =  (1-exp(-u^2/(2*delta^2)))z!phiz(u) = (u^2)/(2*delta^2 + u^2)z"phiz(u) = log(1 + (u^2)/(delta^2))z"phiz(u) =  sqrt(1 + u^2/delta^2)-1zphiz(u) = 1/2 u^2z	lambda = zdelta = zkappa = zeta = zxmin = z and xmax = N)�y�H�H_adj�eta�lambda_�delta�kappa�phi�x�xstar�xbar�xmin�xmax�NbIt�timemax�Crit�Time�NormX�Ndx�SNR�Err�Mem�shape�NxZNyZNz�stop�modaff�print)�selfr   r   r   r   r   r   r   r   r   r   r   r   r    r!   r"   ZphixyZphiz� r0   �G/home/mathieuchalvidal/PycharmProjects/untitled/BP3MG/BP3MG_Standard.py�__init__   sv    																									






z)MajorizeMinimizeMemoryGradient3D.__init__c             C   si  |  j  } |  j } |  j } t | |  j |  j |  j |  j |  j |  j	 |  j
 |  j |  j |  j � \ } } |  j j | � t d | � |  j j t j �  � |  j j t j j | | � � |  j j t j � | d k r� |  j j t j � nX |  j j t j j | | � t j j | � � |  j j d t j d |  j d � � x�t |  j � D]w} |  j j t j �  � |  j | |  j k  r�P| |  j  d k r�t d j! | � � t" | � \ } } }	 t# | |  j d � }
 | d k r�|  j$ | | | | |	 |
 |  j |  j |  j |  j	 |  j
 |  j |  j � } t | j% d � t j j | � | } | | } | |
 } | | } | | } | |	 } n@|  j& | | | | | |	 | | | |
 | |  j |  j |  j |  j	 |  j
 |  j |  j � } t j j | � } t j' | j( d � | j( d � � } t j j) | � t j* | | g � } | d | | d | } | d |
 | d | } | d | | d | } | d | | d | } | d |	 | d | } | | } t | |  j |  j |  j |  j |  j |  j	 |  j
 |  j |  j |  j � \ } } |  j j | � |  j j t j �  � | d k	 r�|  j j t j j | | � t j j | � � |  j j d t j d |  j d � � |  j j t j j | | � � |  j j t j j | � t j j | � � | |  j  d k r�t d	 � t d
 | � t d |  j d � t d |  j d � t d |  j d � t d |  j d |  j d � | d k	 r�t d |  j d � |  j d |  j d |  j+ k rYt d � PqYW|  j d } t d
 t, |  j � � t d |  j d  |  j d � t d |  j d! � t d | � t d � | |  j | |  j |  j |  j |  j- f S)"z&
        Starting computation
        zInitial criterion value = Nr   r   r   u   to be plotted: iteration n°:{}�same�0z---zIteration number = zCriterion value = zError to ground truthzNorm(dx) = zComputation time = zSNR value = zMAXIMAL TIME!!zComputation time (cpu) =zFinal criterion value = zFinal SNR value = z(****************************************�����r5   r5   r5   r5   r5   r5   r5   r5   r5   r5   r5   r5   ).r   r   r   r   r   r   r   r   r   r   r   r   r   r    r#   �appendr.   r$   �timer%   �np�linalg�normr&   �infr'   �Infr(   �log10�ranger!   r,   r-   �formatr   r   �
majorante1r*   �
majorante2�dot�reshape�pinv�arrayr"   �lenr)   )r/   r   r   r   ZCriZGrad�kZVvgZVhgZVtgZHg�B�s�dxZHdxZVvdxZVhdxZVtdx�d1�d2�SNRendr0   r0   r1   �optimizep   s�    			N 0(I\$&
N0( ,

!

z)MajorizeMinimizeMemoryGradient3D.optimizec             C   s�  | d } | d } t  | � \ } } } | d k ry d | d t j t j | d | d � d d | d � } n�| d k r� d | d d | d t j | d | d � d d } n�| d k r� d | d t j | d | d � d } nW| d k rGd | d d t j | d | d � d | d d } n| d k rzt j t j | d | d � � } n� | d k r�d | d d t j | d | d � d | d d t j d t j | d | d � d | d d d � } nU | d k rTd	 } d | | d d t j | d | d � d | d | d } d
 } | d k r�d | d t j | d d | d � } n� | d k r�d | d d | d | d d } nu | d k r�d | d | d } nP | d k r$d | d d | d | d d } n | d k r?t j | � } t j | d � |  j t j | | | | | | � |  j t j | | | � |  j t j | | | k | | k d � } | S)Nr   r   r   r   r   r
   r   r   g      �?g�������?r5   g      �r5   g      �g      �?r5   g      �)	r   r8   �exp�sqrt�diag�sumr   r   r   )r/   r   rK   �Vvd1�Vhd1�Vtd1�Hd1r   r   r   r   r   r   r    �phiXY�phiZ�Vvx�Vhx�Vtx�wXY_Vx�pow�p�wZ_VxrH   r0   r0   r1   r@     s>    

D>.>'yC/))	�z+MajorizeMinimizeMemoryGradient3D.majorante1c       -      C   s{  | d } | d } t  | � \ } } } | d k ry d | d t j t j | d | d � d d | d � } n�| d k r� d | d d | d t j | d | d � d d } n�| d k r� d | d t j | d | d � d } nW| d k rGd | d d t j | d | d � d | d d } n| d k rzt j t j | d | d � � } n� | d k r�d | d d t j | d | d � d | d d t j d t j | d | d � d | d d d � } nU | d k rTd	 } d | | d d t j | d | d � d | d | d } d
 } | d k r�d | d t j | d d | d � } n� | d k r�d | d d | d | d d } nu | d k r�d | d | d } nP | d k r$d | d d | d | d d } n | d k r?t j | � } |  j t j | | | | | | � |  j t j | | | � } |  j t j | | | | | | � |  j t j | | |	 � } |  j t j | | | | | | � |  j t j |	 | |	 � } t j | | k | d d � } t j | | k | d d � }  t j | | k | d d � }! t j | | k | d d � }" | }# |! }$ | |  }% |! |" }& |  d }' |" d }( t j |
 d � | |  j	 t j | |! � }) t j |
 | � | |  j	 t j | |  |! |" � }* t j | d � | |  j	 t j |  |" � }+ t j
 |) |* g |* |+ g g � }, |, S)Nr   r   r   r   r   r
   r   r   g      �?g�������?r5   g      �r5   g      �g      �?r5   g      �)r   r8   rO   rP   rQ   r   rR   r   �wherer   rE   )-r/   r   rK   rL   rS   rT   rU   ZVvd2ZVhd2ZVtd2rV   ZHd2r   r   r   r   r   r   r    rW   rX   rY   rZ   r[   r\   r]   r^   r_   Z	d1tVtWVd1Z	d1tVtWVd2Z	d2tVtWVd2Zd1_minZd2_minZd1_maxZd2_maxZ	d1UNmind1Z	d1UNmaxd1Z	d1UNmind2Z	d1UNmaxd2Z	d2UNmind2Z	d2UNmaxd2ZB11ZB12ZB22rH   r0   r0   r1   rA   @  s^    

D>.>'yC/))	FFF



/7/!z+MajorizeMinimizeMemoryGradient3D.majorante2N)�__name__�
__module__�__qualname__r2   rN   r@   rA   r0   r0   r0   r1   r      s   ^�:r   )�numpyr8   r   r   r   r   r   r7   �scipy.signalr   r   �scipy.io�scipy�io�loadmat�brainrE   �Gr   r0   r0   r0   r1   �<module>   s   "3