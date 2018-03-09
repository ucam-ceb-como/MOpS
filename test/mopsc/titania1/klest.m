function kl_ij = klest(dnnew_i,dnnew_j,w_i,w_j)
    m_i = sum(dnnew_i.*w_i)/sum(w_i);
    m_j = sum(dnnew_j.*w_j)/sum(w_j);
    v_i = sum(((dnnew_i-m_i).^2).*w_i)/(sum(w_i)-1);
    v_j = sum(((dnnew_j-m_j).^2).*w_j)/(sum(w_j)-1);
    
    mu_i = log(m_i/sqrt(1+(v_i/(m_i^2))));
    mu_j = log(m_j/sqrt(1+(v_j/(m_j^2))));
    ss_i = log(1+(v_i/(m_i^2)));
    ss_j = log(1+(v_j/(m_j^2)));
    
    kl_ij = (1/(2*ss_j))*((mu_i-mu_j)^2+ss_i-ss_j)+log(sqrt(ss_j/ss_i));